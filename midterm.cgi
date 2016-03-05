#!/usr/local/bin/python3

import jinja2
import re

import cgi
import cgitb

'''
Strategy:
This parses a genBank output file and a glimmer output file comparing start
and end sites.


Input:
genBankFile - a path to the genBankFile output file
glimmerFile - a path to the glimmer output file
    (the .predict file)

Output:
Returns a dictionary of genBank info (directionality/start/stop/protein name) and glimmer
info (directionality/start/stop/orfName) as well as the relationship (or lack 
thereof) between them

'''
def compareGenBankGlimmer(genBankFile, glimmerFile):
    genBankEntries = []
    newEntry = {}
    
    #get the genBank data
    inputFile = open(genBankFile)
    for line in inputFile:
        line = line.rstrip()
        m = re.search(r"\sCDS\s", line)
        if (m):
            CDS = line.split()[1]
            if (CDS.startswith("join")):
                joinCDS = CDS.strip("join(").strip(")")
                components = joinCDS.split("..")
                newEntry['start'] = components[0]
                newEntry['end'] = components[2]    
                newEntry['direction'] = "forward"
                newEntry['Consensus'] = "none"
                newEntry['partner'] = "none"
            elif (CDS.startswith("complement")):
                m = re.search(r"\((\S+)\.\.(\S+)\)", CDS)
                components = m.group(0).replace("(", "").replace(")", "").split("..")
                newEntry['start'] = components[0].strip(">").strip("<")
                newEntry['end'] = components[1].strip(">").strip("<")
                newEntry['direction'] = "complement"
                newEntry['Consensus'] = "none"
            else:
                m = re.search(r"\((\S+)\.\.(\S+)\)", CDS)
                components = CDS.split("..")
                newEntry['start'] = components[0].strip("(").strip(">").strip("<")
                newEntry['end'] = components[1].strip(")").strip(">").strip("<")
                newEntry['direction'] = "forward"
                newEntry['Consensus'] = "none"
        if (re.search("/protein_id=", line)):
            m = line.rstrip()
            m = m.replace("protein_id=", "").replace("/", "").replace("\"", "").strip()
            newEntry['protein_id'] = m
            genBankEntries.append(newEntry)
            newEntry = {}
    inputFile.close()

    #get the predicted glimmer orf data
    predictedEntries = []
    inputFile = open(glimmerFile)
    for line in inputFile:
        newEntry = {}
        line = line.rstrip()
        if (line.startswith("orf")):
            components = line.split()
            newEntry['orfName'] = components[0]
            newEntry['start'] = components[1]
            newEntry['end'] = components[2]
            newEntry['readingFrame'] = int(components[3])
            if (int(components[3]) > 0):
                newEntry['direction'] = "forward"
            else:
                newEntry['direction'] = "complement"
            newEntry['rawScore'] = float(components[4])
            newEntry['Consensus'] = "none"
            predictedEntries.append(newEntry)
    inputFile.close()

    #tie the relationships together
    for x in predictedEntries:
        for y in genBankEntries:
            if ((x['start'] == y['start']) and (x['end'] == y['end']) and (x['direction'] == y['direction'])):
                x['5_match'] = "Agree"
                x['3_match'] = "Agree"
                x['Consensus'] = "Exact Match"
                y['Consensus'] = "Exact Match"
                x['partner'] = y
                y['partner'] = x
            elif ((x['start'] == y['start']) and (x['end'] != y['end']) and (x['direction'] == y['direction'])):
                x['5_match'] = "Agree"
                x['3_match'] = "Disagree"
                x['Consensus'] = "5' Match"
                y['Consensus'] = "5' Match"
                x['partner'] = y
                y['partner'] = x
            elif ((x['start'] != y['start']) and (x['end'] == y['end']) and (x['direction'] == y['direction'])):
                x['5_match'] = "Disagree"
                x['3_match'] = "Agree"
                x['Consensus'] = "3' Match"
                y['Consensus'] = "3' Match"
                x['partner'] = y
                y['partner'] = x
                
    #now go through and assign where there is no relationship
    for x in genBankEntries:
        if (x['Consensus'] == "none"):
            x['Consensus'] = "No Matches"
            x["5_match"] = "Disagree"
            x["3_match"] = "Disagree"
            x['partner'] = "none"
    
    #now go through and assign where 
    for y in predictedEntries:
        if (y['Consensus'] == "none"):
            y['Consensus'] = "No Matches"
            y["5_match"] = "Disagree"
            y["3_match"] = "Disagree"
            y['partner'] = "none"
    return [genBankEntries, predictedEntries]



'''
Strategy:
This calculates out the comparison statistics

Input:
genBankEntries - a list of dictionary of genBankEntries
predictedEntries - a list of dictionary of predictedEntries

Output:
The comparison statistics

'''
def getComparisonData(genBankEntries, predictedEntries):
    #gene count in reference annotation
    genBankCount = len(genBankEntries)
    #gene count in predicted genes
    predictedCount= len(predictedEntries)
    #count of genes with exact matching coordinates 
    exactCount = (sum(x['Consensus'] == "Exact Match" for x in genBankEntries))
    #count of 5' agreement and 3' disagreement
    prime5Count = (sum(x['Consensus'] == "5' Match" for x in genBankEntries))
    #count of 5' disagreement and 3' agreement
    prime3Count = (sum(x['Consensus'] == "3' Match" for x in genBankEntries))
    #genes predicted with no overlap
    genBankNoOverlapCount = genBankCount - exactCount - prime5Count - prime3Count
    predictedNoOverlapCount = predictedCount - exactCount - prime5Count - prime3Count
    return [genBankCount, predictedCount, exactCount, prime5Count, prime3Count, genBankNoOverlapCount, predictedNoOverlapCount]
    
    
'''
Strategy:
This formats the data for the html


Input:
genBankEntries - a list of dictionary of genBankEntries
predictedEntries - a list of dictionary of predictedEntries

Output:
A list of formatted data
 
'''    

def formatOutput(genBankEntries, predictedEntries):
    fullComparison = []
    for x in genBankEntries:
        newEntry = []
        newEntry.append(x['protein_id'])
        #predicted ORF
        if (x['partner'] != "none"):
            newEntry.append(x['partner']['orfName'])
        else:
            newEntry.append("No Match")   
        #start
        newEntry.append(x['start'])
        #glim 5 start
        if (x['partner'] != "none"):
            newEntry.append(x['partner']['start'])
        else:
            newEntry.append("N/A")   
        #5 match
        if (x['partner'] != "none"):
            newEntry.append(x['partner']['5_match'])
        else:
            newEntry.append("N/A")        
        #end
        newEntry.append(x['end'])
        #glim 3 end
        if (x['partner'] != "none"):
            newEntry.append(x['partner']['end'])
        else:
            newEntry.append("N/A") 
        #3 match
        if (x['partner'] != "none"):
            newEntry.append(x['partner']['3_match'])
        else:
            newEntry.append("N/A")        
        #end
        newEntry.append(x['Consensus'])
        fullComparison.append(newEntry)
    return(fullComparison)


templateLoader = jinja2.FileSystemLoader( searchpath="./templates" )
env = jinja2.Environment(loader=templateLoader)
form = cgi.FieldStorage()
template = env.get_template('index.html')
cgitb.enable()
if form.getvalue('dropdown'):
   glimmerFile = form.getvalue('dropdown') + ".predict"
else:
   glimmerFile = "eColi_o30_g100_t30.predict"
genBankFile = "AB011549.gb"

compare = compareGenBankGlimmer(genBankFile, glimmerFile)
genBankEntries = compare[0]
predictedEntries = compare[1]
data = getComparisonData(genBankEntries, predictedEntries)
fullComparison = formatOutput(genBankEntries, predictedEntries)

genBankCount = data[0]
predictedCount = data[1]
exactCount = data[2]
prime5Count = data[3]
prime3Count = data[4]
genBankNoOverlapCount = data[5]
predictedNoOverlapCount = data[6]

print("Content-Type: text/html\n\n")
print(template.render(values = fullComparison, genBankCount = genBankCount, predictedCount = predictedCount, exactCount = exactCount, prime5Count = prime5Count, prime3Count = prime3Count, genBankNoOverlapCount = genBankNoOverlapCount, predictedNoOverlapCount = predictedNoOverlapCount))








