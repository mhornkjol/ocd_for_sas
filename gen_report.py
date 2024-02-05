import glob 

head = r"""
\documentclass{article}
\usepackage{graphicx,subcaption}
\begin{document}

\textwidth=12cm
"""

tail = """ 
\end{document}
"""

pat_str = """
\\begin{figure}
\\begin{center}
\\begin{subfigure}{0.45\\textwidth}
\includegraphics[width=0.95\\textwidth]{%(file1)s}
\caption{%(file1)s }
\end{subfigure}
\\begin{subfigure}{0.45\\textwidth}
\includegraphics[width=0.95\\textwidth]{%(file0)s}
\caption{%(file0)s }
\end{subfigure}
\\begin{subfigure}{0.45\\textwidth}
\includegraphics[width=0.95\\textwidth]{%(file2)s}
\caption{%(file2)s }
\end{subfigure}
\\begin{subfigure}{0.45\\textwidth}
\includegraphics[width=0.95\\textwidth]{%(file3)s}
\caption{%(file3)s }
\end{subfigure}
\end{center}
\end{figure}
\\newpage
""" 


pats = [105,  175]
pats = [105, 175, 178,  190,  199,  215,  227,  230,  236, 105, 176,  183,  191,  205,  218,  228, 235,  240]

file_str = head 

for pat in pats: 

  pat = f'{pat:03}' 
  fs = glob.glob("%s*png" % pat) 
  fs.sort() 
  d = {}
  for l in range(len(fs)): 
    print (fs[l])
    d["file%d" % l] = fs[l]

  print (d)
  file_str += pat_str % d  
  
file_str += tail 

f = open ("reportDTI.tex", "w") 
f.write(file_str)
f.close()




