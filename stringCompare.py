import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

referenceFile = open('reference.txt','r')
reference = referenceFile.readlines()
referenceFile.close()
for i in range(len(reference)):
    reference[i] = reference[i][:-1]
reference = "".join(reference)

answerFile = open('output.txt','r')
answer = answerFile.readlines()[0]
answerFile.close()
for i in range(len(reference)):
    if reference[i] != answer[i]:
        print("diff at index {}. expetected {} but got {}".format(i,reference[i],answer[i]))