import csv
import argparse


# a class for a tree node with 3 children
class TreeNode(object):
    def __init__(self, left=None, middle=None, right=None):
        self.left = left
        self.middle = middle
        self.right = right


# a function that works recursively and runs the Huffman tree from root to leaves in order to assign
# Huffman codes to the leaves 
def AssigneCodeToChildren(TreeNode, code, huffman_codes):

    l = TreeNode.left
    m = TreeNode.middle
    r = TreeNode.right

    if ( isinstance(l[0], basestring) ):
		x = [l[0], code + "0"]
		huffman_codes.append(x)
    else:
        AssigneCodeToChildren(l[0], code + "0", huffman_codes)
           
    if ( isinstance(m[0], basestring) ):
		x = [m[0], code + "1"]
		huffman_codes.append(x)
    else:
        AssigneCodeToChildren(m[0], code + "1", huffman_codes)
        
    if (( isinstance(r[0], basestring) )):
		x = [r[0], code + "2"]
		huffman_codes.append(x)
    else:
        AssigneCodeToChildren(r[0], code + "2", huffman_codes)
        


# the function that converts the Huffman code to a DNA sequence
def HuffmanCodeToDnaSequence(string):
	
	# incorporate the Huffman-to-DNA transformation table to a special Python dictionary
    huff2dna_keys = ( ("A","0") , ("A","1") , ("A","2") , ("C","0") , ("C","1") , ("C","2") , ("G","0") , ("G","1") , ("G","2") , ("T","0") , ("T","1") , ("T","2") )
    huff2dna_vals = ["C","G","T","G","T","A","T","A","C","A","C","G"]
    huff2dna_dict = dict.fromkeys(huff2dna_keys)
    j = 0
    for i in huff2dna_keys:
		huff2dna_dict[i] = huff2dna_vals[j]
		j = j + 1
	
    huffman_code = string
    dna_seq = "A"
	
	# run the string of numbers and, using the dictionary, match numbers to nucleobases
    for i in range(len(huffman_code)):
        x = dna_seq[i]
        y = huffman_code[i]
        dna_seq = dna_seq + huff2dna_dict.get((x,y))
	
    dna_seq = dna_seq[1:]
    return dna_seq
	


# the function that converts the DNA sequence to a Huffman code
def DnaSequenceToHuffmanCode(string):
	
	# incorporate the DNA-to-Huffman transformation table to a special Python dictionary
    dna2huff_keys = ( ("A","C") , ("A","G") , ("A","T") , ("C","A") , ("C","G") , ("C","T") , ("G","A") , ("G","C") , ("G","T") , ("T","A") , ("T","C") , ("T","G") )
    dna2huff_vals = ["0","1","2","2","0","1","1","2","0","0","1","2"]
    dna2huff_dict = dict.fromkeys(dna2huff_keys)
    j = 0
    for i in dna2huff_keys:
		dna2huff_dict[i] = dna2huff_vals[j]
		j = j + 1
	
	# add an A at the start of the dna sequence
    dna_seq = "A" + string
    huffman_code = ""
    
    # run the dna sequence and, using the dictionary, match nucleobases to numbers 
    for i in range(1, len(dna_seq)):
		x = dna_seq[i-1]
		y = dna_seq[i]
		huffman_code = huffman_code + dna2huff_dict.get((x,y))
		
    return huffman_code


# the encoding function
def encode(input_file, output_file, huffman_file):

    # read input file. Contains raw text
    f = open(input_file, "rb")
    txt = f.read()

    # len_str is the number of characters in the string
    # list_of_chars is a Python list with all the characters in the correct order
    # set_of_chars is a list with the characters that appear at least once
    # pq is the priority queue
    len_str = len(txt)
    list_of_chars = []
    set_of_chars = []   
    pq = []
    
    for i in range(len_str):
        list_of_chars.append(txt[i])

    for i in list_of_chars:
		if i not in set_of_chars:
			set_of_chars.append(i)
    
    # add unique characters with their frequency in the priority queue
    for i in set_of_chars:
        x = [i, list_of_chars.count(i)]
        pq.append(x)
    
    # add a special void character with zero frequency in the priority queue 
    # when the number of different characters in the text is even 
    if ( len(set_of_chars) % 2 == 0 ):
		pq.append(['void',0])

	# sort the priority queue
    pq.sort(key=lambda tup: tup[1])
    
    # create Huffman tree
    while (len(pq) > 1):
        x = pq[0]
        y = pq[1]
        z = pq[2]
        summ = x[1] + y[1] + z[1]
        node = TreeNode(x, y, z)
        pq.insert(0,[node, summ])
        pq.remove(x)
        pq.remove(y)
        pq.remove(z)
        pq.sort(key=lambda tup: tup[1])

    root = pq[0]
    root = root[0]
    root_code = ""
    huffcodes = []

    AssigneCodeToChildren(root, root_code, huffcodes)
    
    # write the Huffman codes to a CSV file
    f = open(huffman_file, "w+")
    huff_writer = csv.writer(f, delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
    huff_writer.writerows(huffcodes)
    f.close()
    
    huff_dict = {}
    for i in huffcodes:
		huff_dict[i[0]] = i[1]

    code_string = ""
    for i in list_of_chars:
		code_string = code_string + huff_dict[i]
	
    dna_seq = HuffmanCodeToDnaSequence(code_string)
    
    f = open(output_file, "w+")
    f.write(dna_seq)
    f.close()
	

# the decoding function
def decode(input_file, output_file, huffman_file):
	
	# read input file. Contains a DNA sequence
    f = open(input_file, "rb")
    dna_seq = f.read()

	# convert dna sequence to a string of numbers
    code_string = DnaSequenceToHuffmanCode(dna_seq)
    
    # read the huffman codes from the huffman file
    f = open(huffman_file, "rb")
    huff_reader = csv.reader(f, delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
    
    # create a dictionary with the huffman codes
    huff_dict = {}
    for row in huff_reader:
		huff_dict[row[1]] = row[0]
    
    # read the string of numbers and assign characters to codes
    # in order to get the original text
    i = 0
    txt = ""
    while (i < len(code_string)):
		j = i + 1
		while (code_string[i:j] not in huff_dict.keys()):
			j = j + 1
		txt = txt + huff_dict[code_string[i:j]]
		i = j

    # save text to the output file
    f = open(output_file, "w+")
    f.write(txt)
    f.close()
		
	
def main():
	
	parser = argparse.ArgumentParser()
	
	# the optional -d argument
	parser.add_argument("-d", "--decode", help="Encode or Decode option.", action="store_true")
	
	# the input file
	parser.add_argument("input", type=str, help="The name of the input file. " +
                        "Should contain text for the encoding option and a DNA sequence" +
                        "for the decoding option")
	
	# the output file
	parser.add_argument("output", type=str, help="The name of the output file. " +
                        "Program outputs a DNA sequence in the encoding option or text " +
                        "in the decoding option")
	
	# the CSV file with the Huffman codes
	parser.add_argument("huffman", type=str, help="The name of the CSV file in which the" +
                        "Huffman coding table is stored")
                      
	args = parser.parse_args()

	if args.decode:
		decode(args.input, args.output, args.huffman)
	else:
		encode(args.input, args.output, args.huffman)


if __name__ == '__main__':
    main()
