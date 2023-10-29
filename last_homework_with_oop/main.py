import sys
from seq2protein import DNA, RNA, Protein, Sequence

'''
Run main.py from command line using
python3 main.py input.fasta output.txt

Input file format:
>sequence_name sequence_type
ATCTCATCTACATC
'''

def main():
    input_file = sys.argv[1] # input.fasta
    output_file = sys.argv[2] # output.txt

    # Your code below, for example
    dna = DNA(input_file)
    res = dna.to_protein()
    with open(output_file, 'w') as out:
        out.write(res)

if __name__ == '__main__':
    main()
