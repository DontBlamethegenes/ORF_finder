# Vanilla Python Gene Finder
By Alberto Aldapa


## Problem: 
Finding all the possible genes in a genome can take a very long time if done incorrectly. Regular expressions fail at edge cases and don't allow for useful user-specified parameters such as minimum length of 
Open Reading Frame ORF.

## Solution:
A python class that takes in a .fasta genome file and returns a text document containing all the potential gene locations. 

## Techniques:
- List of lists of lists
- Class init
- Commandline user parameters 
- Efficient and robust search algorithm
- Printing results to txt file
- Dictionaries
- Tuples

## Results:
![alt text](orfstxt_screenshot.jpg)
Reference | Start Location | Stop Location | Length 

Each potential gene populates one line. 


