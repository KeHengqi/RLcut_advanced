import re
filename = "./output/in-going/wiki/wiki-0.1.csv"
fileopen = open(filename, "r")
file_content = fileopen.readline()
while file_content:
    file_content = fileopen.readline()
    print(file_content)
    print(type(file_content))
    file_split_result = re.split(r'[, \n]', file_content)
    print(file_split_result)
    print(file_split_result[1])
    print(file_split_result[2])
    print(file_split_result[3])
    print(file_split_result[4])
    #print(type(file_split_result))