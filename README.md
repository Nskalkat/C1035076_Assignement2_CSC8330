Advanced Programming for Digital Biology CSC8330
========================================




The purpose of this program:
=====================

The purpose of this code is for a user with a set of DNA sequences in a fastA file to make a consensus sequence and 
then to give the user information based on the consensus sequence. This program will also translate that consensus 
sequence into a protein sequence and from there give information back to user about the protein sequence all in an 
output file. This program is separated into 3 different program files, the consensus_seq.py file, the Protein_Param.py 
file, and the Protein_Param_Main.py file.


How to Operate the Program:
=====================

To operate this code first head to the directory where all the python files are in. Then either on a IDE terminal or 
terminal found on your device type in python Protein_Param_Main.py (has the ability to run all code)the name of your 
FastA file and the consensus threshold value as seen in the example below:

```
python Protein_Param_Main.py test_file_CSC8330.fasta 0.5
```

Now the consensus value can be any value from 0-1 and how the value works is that the closer the consensus value gets 
to 1 the higher your threshold requirement will. Regardless of how high the threshold value is; it will still pick the 
most matched letter in the sequence alignment at that specific position. For example, if there is an alignment of 11 
sequences, in which 6 are A and 5 are C, the A will be chosen because it has a higher match in a sequence. Now when it 
comes to when a match cannot be made, the original program usually gives an ambiguous nucleotide assigned to that 
position, but the problem with that is that an ambiguous will not allow for the rest of this program to analyze this 
sequence and therefore inhibit further use of this program. Since finding the correct nucleotide became a complex 
project by itself, I instead used code to randomly replace an ambiguous N value with a random nucleotide, allowing 
for a sequence that can be analyzed to be created. Once the consensus sequence is taken, it is analyzed for its GC 
content and nucleotide composition. Then the consensus sequence is taken to the Protein_Param.py file where it is 
translated and analyzed using various Protein Analysis techniques and gives back the amino acid composition, 
aromaticity, Isoelectric point etc. This file also takes all the information and returns it to the user as a separate 
file called output.txt that will be created and stored in the same directory as the Protein_Param_Main.py will be 
stored. One caveat of this is that I couldn’t figure out how to return each occurrence as a different 
file, as it will rewrite everything on the same original file. Finally the Protein_Param_Main.py is the python file 
with the Main code that will allow it to run everything using the example code above.

The program also includes an example file called test_file_CSC8330.fasta which contains 10 different DNA sequences and 
can be used to make a consensus sequence as well as give an output file with all the required information

Here is the link to the project on GitHub: 
https://github.com/Nskalkat/C1035076_Assignement2_CSC8330






