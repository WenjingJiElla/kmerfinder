from prettytable import PrettyTable

def readFastaFile(filename):
    """Read the sequence in the Fasta file. For the Fasta file of dsRNA, 
    it should only contain one sequence. For normal sequence Fasta file
    , it could contain more than one sequence. Name cannot be repeated.
    param filename: the fasta file 
    """
    file = open(filename)
    seqlist = {}
    seqName = ""

    for line in file:
        line = line.strip()
        if line[0] == "":
            continue
        elif line[0] == '>':
            seqName = line[1:]
            seqlist[seqName] = ""
        else:
            seqlist[seqName] = line.upper()

    return seqlist

def numKmer(filename, numGene, lenKmer=21):
    """Calculate the number of kmer for dsRNA construct. The length of kmer defaults to 21.
    :param filename: the fasta file of dsRNA construct
    :param numGene: number of genes made up the dsRNA construct
    :param lenKmer: the length of kmer, default is 21 
    """
    seq = list(readFastaFile(filename).items())
    dsRNA = seq[0]

    totalNum = len(dsRNA)-lenKmer+1
    junction = (numGene-1) * (lenKmer-1)
    numKmer = totalNum - junction

    return numKmer

def find_kmers(construct,seq,lenKmer=21):
    """Find the number of kmers that will appear in the target sequence
    :param construct: the sequence of RNA construct
    :param seq: the sequence
    :param kmer_len: the length of kmer, default is 21
    """
    cons = construct.upper() 
    seq = seq.upper() 

    kmers = set() # find all kmers 
    for i in range(len(cons)-lenKmer+1):
        kmers.add(cons[i:i+lenKmer])
     
    count = 0
    for j in kmers:
        if j in seq:
            count +=1
    
    return count

def get_table_col_row(record):
    """Build the table for output.
    :param record: the output of kmer calculation
    """
    # change dict to list
    data = record.items()
    
    list_data = []
    for i in range(len(data)):
        list_data.append(list(data[i])) #[['01000 F_C', 174],...]
    
    # split all item   
    new_list_data = []
    for i in range(len(list_data)):
        new_list_data.append(' '.join(map(str, list_data[i])).split(' ')) #[['01000', 'F_C', '174'],...]

    # create table col and row   
    original_list = []
    target_list = []
    for i in range(len(new_list_data)):
        original_list.append(new_list_data[i][0])
        target_list.append(new_list_data[i][1])
        
    colLabels = list(set(original_list))
    rowLabels = list(set(target_list))
    
    return colLabels, rowLabels

def print_table(record):
    """Print table.
    :param record: the output of kmer calculation
    """
    # original table col and row
    colLabels, rowLabels = get_table_col_row(record)
    colLabels.sort()
    rowLabels.sort()
    
    # get table value
    val = []
    for i in range(len(rowLabels)):
        row = []
        row.append(rowLabels[i])
        for j in range(len(colLabels)):
            index = colLabels[j] + ' '+ rowLabels[i]
            if index in record.keys():
                row.append(record[index])
            else:
                row.append('Null')
        val.append(row)
    # table col
    table_colLabels = colLabels.copy()
    table_colLabels.insert(0, 'SP/G')
    # create a table
    table = PrettyTable(table_colLabels)
    # add table row
    for i in range(len(val)):
        table.add_row(val[i])
    
    print(table)

def matchedKmer(dsRNA_file, sequences_file, numGene):
    """Find the number of kmer appear in the different species.
    :param dsRNA_file: the sequence of RNA construct
    :param sequences_file: the sequence
    :param numGene: the length of kmer, default is 21
    """
   #obtain the info of dsRNA construct
    dsRNA = readFastaFile(dsRNA_file)
    if len(dsRNA) != 1:
        raise Exception("Wrong dsRNA")
    else:
        consName = list(dsRNA.keys())[0]
        consSeq = list(dsRNA.items())[0][1]
    
    print("The total(or Maximum) Kmar for ", consName, "should be:", numKmer(consSeq, numGene)) #total num of kmer  
    
    # obtain the info of sequences
    sequences = readFastaFile(sequences_file)
    # numSpecies = len(sequences) # number of kmer appear in each sequence

    # get the matched kmar for each sequences
    record = sequences.copy() 
    for i in record:
        record[i] = find_kmers(consSeq,sequences[i]) 
    
    return print_table(record)

    





    
        



