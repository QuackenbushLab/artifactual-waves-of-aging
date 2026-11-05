import sys
filename = sys.argv[1]
file = open(filename,"r")
text = file.read().split("\n")
bigstring = "::".join(text)
entries = bigstring.split("::::")
for entry in entries:
    entry = entry.split("::")
    title = entry[1].rstrip(",")
    year = entry[4].rstrip(",")
    doi = entry[7].rstrip(",").lstrip("https://doi.org/").rstrip(".")
    print("\t".join([title,year,doi]))
