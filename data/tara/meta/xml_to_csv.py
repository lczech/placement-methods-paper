#!/usr/bin/python

from xml.etree import ElementTree
import os

fields = []
table = []

def uniq(input):
    output = []
    for x in input:
        if x not in output:
            output.append(x)
    return output

for filename in os.listdir('xml'):

    attribs = {}
    attribs["File"] = filename
    fields.append("File")

    with open('xml/' + filename, 'rt') as f:
        tree = ElementTree.parse(f)

        attribs["Accession"] = tree.find("SAMPLE").attrib["accession"]
        attribs["Alias"] = tree.find("SAMPLE").attrib["alias"]
        fields.append("Accession")
        fields.append("Alias")

        # for node in tree.iter():
        #     print node.tag, node.attrib

        for node in tree.findall('.//SAMPLE_ATTRIBUTE'):
            tag = ""
            val = ""
            uni = ""
            for dat in node.iter():
                # print dat.tag, dat.attrib, dat.text
                if dat.tag == "TAG":
                    tag = dat.text
                if dat.tag == "VALUE":
                    val = dat.text.encode('utf-8')
                if dat.tag == "UNITS":
                    uni = dat.text.encode('utf-8')
            attribs[tag] = val
            fields.append(tag)
            if not uni == "":
                attribs[tag+" Unit"] = uni
                fields.append(tag+" Unit")

    table.append( attribs )
    # print attribs
    # for item in fields:
    #     print item, attribs[item]

fields = uniq(fields)

with open('data.csv', 'wt') as f:
    for item in fields:
        f.write(item + "\t")
    f.write("\n")
    for line in table:
        for item in fields:
            # print item, line[item]
            f.write(line[item] + "\t")
        f.write("\n")
