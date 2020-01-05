import re
import sys

fileReader = open(str(sys.argv[1]),'r')

words = fileReader.read()

print(words)


#print(fileReader.read())

#print("Working arrays")
#arrays = re.findall('\[.*\]',fileReader.read())

#print(arrays)

#line = fileReader.readline()

"""
while(line):
    print(line)
    line = fileReader.readline()
"""
