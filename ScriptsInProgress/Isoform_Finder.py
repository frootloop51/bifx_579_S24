## This program reads a fasta file and obtains the unique sequences, the number
## of each unique sequence and produces a new fasta file with only the unique 
## sequences that are annotated with the number of times that sequence occurred
################# Declare Functions #################
# check type of sequence, returns True is protein sequence #
def checkType(clip):            
    ntch = 'actgnACTGN'
    aach = 'defhiklmpqrsvwxyDEFHIKLMPQRSVWXY'
    for l in clip:
        if l in ntch:
            return False
        elif l in aach:
            return True

# get a list of codons from the sequence #
def get_codons(gene):
    codons = []
    start, stop, step = 0, 3, 1
    while stop <= len(gene):
        codon = gene[slice(start, stop, step)]
        codons += [codon]
        start += 3
        stop += 3
    return codons

# translate codons into a peptide #
def translate(codons):
    ucode = {
    'ttt':'F', 'ttc':'F', 'tta':'L', 'ttg':'L', 'tct':'S',
    'tcc':'S', 'tca':'S', 'tcg':'S', 'tat':'Y', 'tac':'Y',
    'tgt':'C', 'tgc':'C', 'tgg':'W', 'ctt':'L', 'ctc':'L',
    'cta':'L', 'ctg':'L', 'cct':'P', 'ccc':'P', 'cca':'P',
    'ccg':'P', 'cat':'H', 'cac':'H', 'caa':'Q', 'cag':'Q',
    'cgt':'R', 'cgc':'R', 'cga':'R', 'cgg':'R', 'att':'I',
    'atc':'I', 'ata':'I', 'atg':'M', 'act':'T', 'acc':'T',
    'aca':'T', 'acg':'T', 'aat':'N', 'aac':'N', 'aaa':'K',
    'aag':'K', 'agt':'S', 'agc':'S', 'aga':'R', 'agg':'R',
    'gtt':'V', 'gtc':'V', 'gta':'V', 'gtg':'V', 'gct':'A',
    'gcc':'A', 'gca':'A', 'gcg':'A', 'gat':'D', 'gac':'D',
    'gaa':'E', 'gag':'E', 'ggt':'G', 'ggc':'G', 'gga':'G',
    'ggg':'G', 'taa':'X', 'tag':'X', 'tga':'X'}
    pep_seq = []
    for codon in codons:
        if codon in ucode.keys():
            pep_seq += [ucode[codon]]
        else:
            pep_seq += 'o'
            
    peptides = "".join([str(i.split()) for j,i in enumerate(pep_seq)])
    return peptides


################# Main Program #################
count = 0
temp = {}
isoforms = {}
reads = 'sample_data.nt.txt'                                        #read file
with open(reads, 'r') as seq, open('out.txt', 'w') as out:
    lines = [line.rstrip() for line in seq]
    for n, line in enumerate(lines):                #determine type of sequence
        if not (line.find('>')==0) and checkType(line) == False:
            trans = translate(get_codons(line))
            trans = (trans.translate({ord(i): None for i in '[]\''}))
            lines[n] = trans
    
    count = {i:lines.count(i) for i in lines[1::2]}    #create count dictionary
            
    for u in range(0,len(lines)-1,2):               #create sequence dictionary
        temp[lines[u]] = lines[u+1:u+2]              
               
    for key, value in temp.items():                      #find unique sequences
        if value not in isoforms.values():
            isoforms[key[-1:]] = value
            
    for n,key,value in zip(count,isoforms.keys(),isoforms.values()):    #output
        num = count[n]
        out.write('>'+str(key)+'  ['+str(num)+']  '+'\n'+str(value)[2:-2]+'\n')