import regex
import sys

fileReader = open(str(sys.argv[1]),'r')

line = fileReader.readline()

while(line):
    print(line)
    line = fileReader.readline()
