import pysam
import HTSeq
import re
import math
import os
import collections
import argparse

def calc_doench_score(seq):
    '''
    Calculate on-target activity score of sgRNA by Cas9 genome editing
    Input seq: 30-mer (4bp + 20bp sgRNA + 3bp PAM + 3bp)
    Output: score
    
    Score and parameters were obtained from original papter as follows:
     - https://media.nature.com/original/nature-assets/nbt/journal/v32/n12/extref/nbt.3026-S1.pdf
     - https://github.com/maximilianh/crisporWebsite/blob/master/doenchScore.py
    '''
    
    params = [
        (1,'G',-0.2753771),   (2,'A',-0.3238875),   (2,'C',0.17212887),   (3,'C',-0.1006662),
        (4,'C',-0.2018029),   (4,'G',0.24595663),   (5,'A',0.03644004),   (5,'C',0.09837684),
        (6,'C',-0.7411813),   (6,'G',-0.3932644),   (11,'A',-0.466099),   (14,'A',0.08537695),
        (14,'C',-0.013814),   (15,'A',0.27262051),  (15,'C',-0.1190226),  (15,'T',-0.2859442),
        (16,'A',0.09745459),  (16,'G',-0.1755462),  (17,'C',-0.3457955),  (17,'G',-0.6780964),
        (18,'A',0.22508903),  (18,'C',-0.5077941),  (19,'G',-0.4173736),  (19,'T',-0.054307),
        (20,'G',0.37989937),  (20,'T',-0.0907126),  (21,'C',0.05782332),  (21,'T',-0.5305673),
        (22,'T',-0.8770074),  (23,'C',-0.8762358),  (23,'G',0.27891626),  (23,'T',-0.4031022),
        (24,'A',-0.0773007),  (24,'C',0.28793562),  (24,'T',-0.2216372),  (27,'G',-0.6890167),
        (27,'T',0.11787758),  (28,'C',-0.1604453),  (29,'G',0.38634258),  (1,'GT',-0.6257787),
        (4,'GC',0.30004332),  (5,'AA',-0.8348362),  (5,'TA',0.76062777),  (6,'GG',-0.4908167),
        (11,'GG',-1.5169074), (11,'TA',0.7092612),  (11,'TC',0.49629861), (11,'TT',-0.5868739),
        (12,'GG',-0.3345637), (13,'GA',0.76384993), (13,'GC',-0.5370252), (16,'TG',-0.7981461),
        (18,'GG',-0.6668087), (18,'TC',0.35318325), (19,'CC',0.74807209), (19,'TG',-0.3672668),
        (20,'AC',0.56820913), (20,'CG',0.32907207), (20,'GA',-0.8364568), (20,'GG',-0.7822076),
        (21,'TC',-1.029693),  (22,'CG',0.85619782), (22,'CT',-0.4632077), (23,'AA',-0.5794924),
        (23,'AG',0.64907554), (24,'AG',-0.0773007), (24,'CG',0.28793562), (24,'TG',-0.2216372),
    ]
    
    intercept =  0.59763615
    gc_high = -0.1665878
    gc_low = -0.2026259

    score = intercept
    
    guide_seq = seq[4:24]
    gc_count = guide_seq.count("G") + guide_seq.count("C")
    
    if gc_count <= 10:
        gc_weight = gc_low
    if gc_count > 10:
        gc_weight = gc_high
        
    score += abs(10-gc_count)*gc_weight
    
    for pos, model_seq, weight in params:
        sub_seq = seq[pos:pos+len(model_seq)]
        if sub_seq == model_seq:
            score += weight
    return 1.0/(1.0+math.exp(-score))

def find_sgRNA_in_polyc_regoin(fasta, db):
    '''
    Search polyC region that can be targeted by spCas9 (PAM is NGG)
    '''
    
    p = re.compile(r'C{6}[ATGC]{14}[ATGC][G]{2}')
    result = collections.namedtuple('PolycGuideRnaResult', ['chr', 'start', 'end', 'guide', 'PAM', 'score', 'is_exon'])
    with pysam.FastxFile(fasta) as fh:
        for entry in fh:
            for m in p.finditer(entry.sequence):
                start = m.start()
                end = m.end()
                score_seq = entry.sequence[start-4:end+3]
                score = calc_doench_score(score_seq)
                seed_seq = entry.sequence[start+6:end-3]
                sgRNA = entry.sequence[start:end]
                pam = sgRNA[-3:]
                if filter_homopolymer(seed_seq):
                    query_iv = HTSeq.GenomicInterval(entry.name, start, end, '+')
                    is_exon_overlapped = find_exon(query_iv, db)
                    yield result(entry.name, start, end, sgRNA, pam, score, is_exon_overlapped)

def filter_homopolymer(seq):
    patterns = ['AAA', 'GGG', 'TTT', 'CCC']
    for p in patterns:
        if p in seq:
            return False # homo-polymer
    else:
        return True

def generate_gff_db(path):
    db = {}
    for gff in HTSeq.GFF_Reader(path):
        if gff.type == 'exon' and gff.iv.strand == '+':
            db[gff.attr['transcript_id']] = gff
    return db

def find_exon(q_iv, db):
    '''
    To test whether sgRNA candidate is overlapped with known exon of genes
    '''
    
    for k, v in db.items():
        if v.iv.overlaps(q_iv):
            return k
    else:
        return False
        

def main(args):
    gff_db = generate_gff_db(args.gff)
    for fasta in os.listdir(args.fasta_dir):
        for sg in find_sgRNA_in_polyc_regoin(os.path.join(args.fasta_dir, fasta), gff_db):
            if sg.score > args.cutoff:
                print '%s,%d,%d,%s,%s,%f,%s' % (sg.chr, sg.start, sg.end, sg.guide, sg.PAM, sg.score, sg.is_exon)

                
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find sgRNA candidates in the polyC region')
    parser.add_argument('--fasta_dir', action='store', required=True, type=str, help='Fasta directly')
    parser.add_argument('--cutoff', action='store', default=0.1, type=float, help='Min threshold of sgRNA activity score [Default: 0.1]')
    parser.add_argument('--gff', action='store', required=True, default=0.1, type=str, help='Path to GFF file')
    args = parser.parse_args()
    main(args)
    
