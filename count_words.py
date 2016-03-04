from collections import Counter
with open("/home/astretton/Downloads/chemistry.dic") as chem:
    words = [line.strip() for line in chem.readlines() if line.strip()]

prefixes = Counter()
suffixes = Counter()

for index, line in enumerate(words):

    for line2 in words:
        if line2:
            if line.endswith(line2) and line != line2:
                suffixes.update([line2])
                prefixes.update([line[0:(len(line)-len(line2))]])
    if index % 100 == 0:
        print "%d/%d" % (index, len(words)) 

from pprint import pprint
pprint(prefixes.most_common(100))
pprint(suffixes.most_common(100))
