Last login: Fri Mar  6 09:52:34 on ttys001
(base) mirna@MacBookPro ~ % git --version
git version 2.50.1 (Apple Git-155)
(base) mirna@MacBookPro ~ % cd /Users/mirna/Desktop/Mirna/MasteOmicsDataAnalysis_Study/Epigenomics/git_HandsOn
(base) mirna@MacBookPro git_HandsOn % git init
Initialized empty Git repository in /Users/mirna/Desktop/Mirna/MasteOmicsDataAnalysis_Study/Epigenomics/git_HandsOn/.git/
(base) mirna@MacBookPro git_HandsOn % git status
On branch main

No commits yet

nothing to commit (create/copy files and use "git add" to track)
(base) mirna@MacBookPro git_HandsOn % git add <file>
zsh: parse error near `\n'
(base) mirna@MacBookPro git_HandsOn % cd/           
zsh: no such file or directory: cd/
(base) mirna@MacBookPro git_HandsOn % cd \
> 
(base) mirna@MacBookPro ~ % cd /Users/mirna/Desktop/Mirna/MasteOmicsDataAnalysis_Study/Epigenomics/git_HandsOn         
(base) mirna@MacBookPro git_HandsOn % nano









  UW PICO 5.09                                           New Buffer                                           Modified  

#!/usr/bin/env python

import sys, re
from argparse import ArgumentParser

parser = ArgumentParser(description = 'Classify a sequence as DNA or RNA')
parser.add_argument("-s", "--seq", type = str, required = True, help = "Input sequence")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

if re.search('^[ACGTU]+$', args.seq):
    if re.search('T', args.seq):
        print ('The sequence is DNA')
    elif re.search('U', args.seq):
        print ('The sequence is RNA')
    else:
        print ('The sequence can be DNA or RNA')
else:
    print ('The sequence is not DNA nor RNA')



^G Get Help         ^O WriteOut         ^R Read File        ^Y Prev Pg          ^K Cut Text         ^C Cur Pos          
^X Exit             ^J Justify          ^W Where is         ^V Next Pg          ^U UnCut Text       ^T To Spell        
