import tkinter as tk
from tkinter import messagebox
from tkinter import filedialog
from Bio import AlignIO
from Bio import SeqIO
from Bio.pairwise2 import format_alignment
from pathlib import Path
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
import os
from Bio import Phylo
import glob
import screed
from Bio import pairwise2
from Bio.SeqUtils import GC
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Entrez

root = tk.Tk()
root.geometry('1500x700')
root.title('SWE')
root.resizable(0,0)


def deletefile():
    
    dir = "/Users/Aya Abo Bakr/Downloads/"

    for path in glob.iglob(os.path.join(dir, '*.aln')):
        os.remove(path)
    for path in glob.iglob(os.path.join(dir, '*.dnd')):
        os.remove(path)
    for path in glob.iglob(os.path.join(dir, '*.phy')):
        os.remove(path)


def seq_MSA():
    deletefile()
    seq = output.get("1.0", "end-1c")
    oh = open(r"C:\Users\Aya Abo Bakr\Downloads\Seq1.fasta",'w')
    oh.write(seq)
    oh.close()
    do_MSA(r"C:\Users\Aya Abo Bakr\Downloads\Seq1.fasta")

def seq_tree():
    deletefile()
    seq = outputtree.get("1.0", "end-1c")
    oh = open(r"C:\Users\Aya Abo Bakr\Downloads\Seq1.fasta",'w')
    oh.write(seq)
    oh.close()
    global phylo
    phylo= r"C:\Users\Aya Abo Bakr\Downloads\Seq1.fasta"


def html():
    deletefile()

    path=filedialog.askopenfilename()
    output.delete('1.0', tk.END)
    cline = MuscleCommandline(input=path, out=r"C:\Users\Aya Abo Bakr\Downloads\muscle.aln", html=True)
    print(cline)
    muscle_exe = r"C:\Users\Aya Abo Bakr\Downloads\muscle3.8.31_i86win32.exe"
    muscle_cline = MuscleCommandline(muscle_exe, input=path, out=r"C:\Users\Aya Abo Bakr\Downloads\muscle.aln", html=True)
    assert os.path.isfile(muscle_exe), "Muscle executable missing"
    stdout, stderr = muscle_cline()
    os.startfile(os.path.abspath(r"C:\Users\Aya Abo Bakr\Downloads\muscle.aln")) 




def do_ncbi(database):
    Entrez.email = "ayaabobakr.uni2002@gmail.com"  
    seqID = seqoutput.get("1.0", "end-1c")
    if(database=="nucleotide"):
        name = seqID+"n"
    if(database=="protein"):
        name = seqID+"p"

    filename = f"{name}.fasta"
    fullpath = os.path.join(r"C:\Users\Aya Abo Bakr\Downloads",filename)
    
    if not os.path.isfile(fullpath):
        net_handle = Entrez.efetch(
          db=database, id=[seqID], rettype="fasta", retmode="text")
        out_handle = open(fullpath, "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        for record in SeqIO.parse(fullpath, "fasta"):
            idtext.delete('1.0', tk.END)
            nametext.delete('1.0', tk.END)
            descriptionttext.delete('1.0', tk.END)
            seqtext.delete('1.0', tk.END)
            idtext.insert(tk.INSERT,record.id)
            nametext.insert(tk.INSERT,record.name)
            seqtext.insert(tk.INSERT,record.seq)
            descriptionttext.insert(tk.INSERT,record.description)
    else:
        messagebox.showinfo("Error", "Already Exists!")
        seqoutput.delete('1.0', tk.END)
    

def do_blast():
    firstseqoutput.config(state="normal")
    secondseqoutput.config(state="normal")
    pwoutput.config(state="normal")
    alignments = pairwise2.align.globalxx(firstseqoutput.get("1.0", "end-1c"),secondseqoutput.get("1.0", "end-1c"))
    pwoutput.insert(tk.INSERT,format_alignment(*alignments[0]))
    firstseqoutput.config(state="disable")
    secondseqoutput.config(state="disable")
    pwoutput.config(state="disable")

    



def do_MSA(file):
    output.config(state="normal")
    output.delete('1.0', tk.END)
    consereved.delete('1.0', tk.END)
    weak_sim.delete('1.0', tk.END)
    strong_sim.delete('1.0', tk.END)
    gaps.delete('1.0', tk.END)
    deletefile()
    cline = ClustalwCommandline("clustalw2", infile=file )
    print(cline)
    clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
    clustalw_cline = ClustalwCommandline(clustalw_exe, infile=file)
    assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
    stdout, stderr = clustalw_cline()
    sp = Path(file).stem
    print(sp)
    full_path=os.path.join("C:\\Users\\Aya Abo Bakr\\Downloads\\",f"{sp}.aln")
    print(full_path)
    if os.path.isfile(full_path):
        oh = open(full_path)
        s = oh.read()[42:]
        output.insert(tk.INSERT, s)
        output.config(state="disabled")
        gaps.insert(tk.INSERT,s.count("-"))
        weak_sim.insert(tk.INSERT,s.count("."))
        strong_sim.insert(tk.INSERT,s.count(":"))
        consereved.insert(tk.INSERT,s.count("*"))
    if os.path.isfile(r"C:\Users\Aya Abo Bakr\Downloads\Seq1.aln"):
        oh = open(r"C:\Users\Aya Abo Bakr\Downloads\Seq1.aln")
        s = oh.read()[42:]
        output.insert(tk.INSERT, s)
        output.config(state="disabled")
        gaps.insert(tk.INSERT,s.count("-"))
        weak_sim.insert(tk.INSERT,s.count("."))
        strong_sim.insert(tk.INSERT,s.count(":"))
        consereved.insert(tk.INSERT,s.count("*"))
    

def do_dnd():
    deletefile()
    outputtree.config(state="normal")
    cline = ClustalwCommandline("clustalw2", infile=phylo )
    print(cline)
    clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
    clustalw_cline = ClustalwCommandline(clustalw_exe, infile=phylo)
    assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
    stdout, stderr = clustalw_cline()
    sp = Path(phylo).stem
    print(sp)
    full_path=os.path.join("C:\\Users\\Aya Abo Bakr\\Downloads\\",f"{sp}.dnd")
    print(full_path)
    if os.path.isfile(full_path):
        tree = Phylo.read(full_path, "newick")
        Phylo.draw(tree)
        outputtree.config(state="disabled")
    if os.path.isfile(r"C:\Users\Aya Abo Bakr\Downloads\Seq1.dnd"):
        tree = Phylo.read(r"C:\Users\Aya Abo Bakr\Downloads\Seq1.dnd", "newick")
        Phylo.draw(tree)
        outputtree.config(state="disabled")

def do_phy():
    deletefile()
    outputtree.config(state="normal")
    cline = MuscleCommandline(input=phylo, out=r"C:\Users\Aya Abo Bakr\Downloads\muscle.phy", phyi=True)
    print(cline)
    muscle_exe = r"C:\Users\Aya Abo Bakr\Downloads\muscle3.8.31_i86win32.exe"
    muscle_cline = MuscleCommandline(muscle_exe, input=phylo, out=r"C:\Users\Aya Abo Bakr\Downloads\muscle.phy",phyi=True)
    assert os.path.isfile(muscle_exe), "Muscle executable missing"
    stdout, stderr = muscle_cline()
    mphy = AlignIO.read(r"C:\Users\Aya Abo Bakr\Downloads\muscle.phy", 'phylip')
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(mphy)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    Phylo.draw(tree)
    outputtree.config(state="disabled")


def browseFiles(n):
    deletefile()
    if(n==3):
        outputtree.delete('1.0', tk.END)
        global phylo
        phylo = filedialog.askopenfilename()
        return
    path=filedialog.askopenfilename()
    if(n==1):
        do_MSA(path)
    if(n==2):
        do_Stat(path)
    



def seq_stat():
    seq = seqoutput.get("1.0", "end-1c")
    oh = open(r"C:\Users\Aya Abo Bakr\Downloads\Seq1.fasta",'w')
    oh.write(seq)
    oh.close()
    do_Stat(r"C:\Users\Aya Abo Bakr\Downloads\Seq1.fasta")


def readFastaFile(inputfile):
    """
    Reads and returns file as FASTA format with special characters removed.
    """
    with screed.open(inputfile) as seqfile:
        for read in seqfile:
            seq = read.sequence
    return seq

def do_Stat(file):
    seqoutput.config(state="normal")
    seqoutput.delete('1.0', tk.END)
    gc.delete('1.0', tk.END)
    bases.delete('1.0', tk.END)
    length.delete('1.0', tk.END)
    s = readFastaFile(file)
    seqoutput.insert(tk.INSERT, s)
    seqoutput.config(state="disabled")
    length.insert(tk.INSERT,len(s))
    base_counts = {} 
    for base in s: 
        if base in base_counts: 
            base_counts[base] += 1 
        else: 
            base_counts[base] = 1 
    bases.insert(tk.INSERT,base_counts)
    gc.insert(tk.INSERT,GC(s))
    
    


def delete_page():
    for frame in main_frame.winfo_children():
        frame.destroy()


def hide_ind():
    ncbi_indicate.config(bg='#1E3694')
    blast_indicate.config(bg='#1E3694')
    msa_indicate.config(bg='#1E3694')
    phy_indicate.config(bg='#1E3694')
    stat_indicate.config(bg='#1E3694')


def indicate(lb, page):
    hide_ind()
    lb.config(bg='#91A1CD')
    delete_page()
    page()



frame = tk.Frame(root, bg='#1E3694')

#NCBI

ncbi_btn = tk.Button(frame, text='NCBI', font=('Bold', 20),
                        fg='White', bd=0, bg='#91A1CD', command=lambda:indicate(ncbi_indicate, ncbi))

ncbi_btn.place(x=60,y=100)

ncbi_indicate = tk.Label(frame, text='', bg='#1E3694')
ncbi_indicate.place(x=3,y=100, width=5, height=55)

def ncbi():

    accession = tk.Label(main_frame,text="Enter the accession number here",bg='#5B72B4', fg = 'Black', font=('Bold', 20))
    accession.place(x=100,y=50)
    Nucleotide_btn = tk.Button(main_frame, text='Nucleotide', font=('Bold', 20),
                        fg='White', bd=0, bg='#1E3694', command=lambda:do_ncbi("nucleotide"))
    Nucleotide_btn.place(x=1000,y=100)

    Protein_btn = tk.Button(main_frame, text='Protein', font=('Bold', 20),
                        fg='White', bd=0, bg='#1E3694', command=lambda:do_ncbi("protein"))
    Protein_btn.place(x=1000,y=200)

    global seqoutput
    seqoutput = tk.Text(main_frame, height = 1,
              width = 50,
              bg = "#91A1CD",font=('Bold', 19))
        
    seqoutput.place(x=100,y=100)

    idlabel = tk.Label(main_frame,text="~ID:",bg='#5B72B4', fg = 'Black', font=('Bold', 15))
    idlabel.place(x=100,y=400)

    namelabel = tk.Label(main_frame,text="~Name:",bg='#5B72B4', fg = 'Black', font=('Bold', 15))
    namelabel.place(x=100,y=500)

    desclabel = tk.Label(main_frame,text="~Description:",bg='#5B72B4', fg = 'Black', font=('Bold', 15))
    desclabel.place(x=600,y=400)

    seqlabel = tk.Label(main_frame,text="~Sequence:",bg='#5B72B4', fg = 'Black', font=('Bold', 15))
    seqlabel.place(x=600,y=500)
 
    global idtext,nametext,descriptionttext,seqtext

    idtext = tk.Text(main_frame, height = 2,
              width = 20,
              bg = "#91A1CD")
        
    idtext.place(x=200,y=400)

    nametext = tk.Text(main_frame, height = 2,
              width = 20,
              bg = "#91A1CD")
        
    nametext.place(x=200,y=500)

    descriptionttext = tk.Text(main_frame, height = 2,
              width = 40,
              bg = "#91A1CD")
        
    descriptionttext.place(x=750,y=400)

    seqtext = tk.Text(main_frame, height = 5,
              width = 40,
              bg = "#91A1CD")
        
    seqtext.place(x=750,y=500)







#Blast
blast_btn = tk.Button(frame, text='BLAST', font=('Bold', 20),
                        fg='White', bd=0, bg='#91A1CD', command=lambda:indicate(blast_indicate, blast))

blast_btn.place(x=50,y=200)

blast_indicate = tk.Label(frame, text='', bg='#1E3694')
blast_indicate.place(x=3,y=200, width=5, height=55)

def blast():

    global pwoutput, secondseqoutput,firstseqoutput

    firstseq = tk.Label(main_frame,text="First Sequence!",bg='#5B72B4', fg = 'White', font=('Bold', 20))
    firstseq.place(x=500,y=50)

    firstseqoutput = tk.Text(main_frame, height = 4,
              width = 125,
              bg = "#91A1CD")
        
    firstseqoutput.place(x=20,y=100)

    secondseq = tk.Label(main_frame,text="Second Sequence!",bg='#5B72B4', fg = 'White', font=('Bold', 20))
    secondseq.place(x=500,y=250)

    secondseqoutput = tk.Text(main_frame, height = 4,
              width = 125,
              bg = "#91A1CD")
        
    secondseqoutput.place(x=20,y=300)
    
    blastbtn = tk.Button(main_frame, text='BLAST!', font=('Bold', 20),
                        fg='White', bd=0, bg='#1E3694', command=lambda:do_blast())
    blastbtn.place(x=550,y=400)

    pwoutput = tk.Text(main_frame, height = 10,
              width = 125,
              bg = "#91A1CD")
        
    pwoutput.place(x=20,y=470)





#MSA

msa_btn = tk.Button(frame, text='MSA', font=('Bold', 20),
                        fg='White', bd=0, bg='#91A1CD', command=lambda:indicate(msa_indicate, msa))

msa_btn.place(x=60,y=300)

msa_indicate = tk.Label(frame, text='', bg='#1E3694')
msa_indicate.place(x=3,y=300, width=5, height=55)

def msa():

    choose_file = tk.Label(main_frame,text="Export a fasta file",bg='#5B72B4', fg = 'Black', font=('Bold', 15))
    choose_file.place(x=10,y=30)

    writeseq = tk.Label(main_frame,text="or write your sequences here",bg='#5B72B4', fg = 'Black', font=('Bold', 10))
    writeseq.place(x=10,y=75)

    align = tk.Button(main_frame,text='Align!',command=lambda:seq_MSA(), font=('Bold', 10),
                        fg='White', bd=0, bg='#91A1CD')
    align.place(x=260,y=70)

    file_btn = tk.Button(main_frame, text='Browse', font=('Bold', 10), width=20,
                        fg='White', bd=0, bg='#91A1CD', command=lambda:browseFiles(1))
    file_btn.place(x=200,y=35)

    html_btn = tk.Button(main_frame,text="Click here for HTML view",fg='White', bd=0, bg='#91A1CD', 
                        font=('Bold', 15), command=lambda:html())
    html_btn.place(x=330,y=650)

    conserevedlbl = tk.Label(main_frame,text="~Conserved Residues:",bg='#5B72B4', fg = 'Black', font=('Bold', 15))
    conserevedlbl.place(x=930,y=100)

    gapslbl = tk.Label(main_frame,text="~Gaps:",bg='#5B72B4', fg = 'Black', font=('Bold', 15))
    gapslbl.place(x=930,y=200)

    weak_simlbl = tk.Label(main_frame,text="~Weakily Similar Residues:",bg='#5B72B4', fg = 'Black', font=('Bold', 15))
    weak_simlbl.place(x=930,y=300)

    strong_simlbl = tk.Label(main_frame,text="~Strongly Similar Residues:",bg='#5B72B4', fg = 'Black', font=('Bold', 15))
    strong_simlbl.place(x=930,y=400)

    global output 
    global weak_sim
    global strong_sim
    global gaps
    global consereved

    consereved = tk.Text(main_frame, height = 2,
              width = 4,
              bg = "#91A1CD")
        
    consereved.place(x=1200,y=90)

    gaps = tk.Text(main_frame, height = 2,
              width = 4,
              bg = "#91A1CD")
        
    gaps.place(x=1200,y=190)

    strong_sim = tk.Text(main_frame, height = 2,
              width = 4,
              bg = "#91A1CD")
        
    strong_sim.place(x=1200,y=290)

    weak_sim = tk.Text(main_frame, height = 2,
              width = 4,
              bg = "#91A1CD")
        
    weak_sim.place(x=1200,y=390)
 
    output = tk.Text(main_frame, height = 30,
              width = 90,
              bg = "#91A1CD")
        
    output.place(x=10,y=100)





#Phylo

phy_btn = tk.Button(frame, text='Phylogenetic\nTree', font=('Bold', 20),
                        fg='White', bd=0, bg='#91A1CD', command=lambda:indicate(phy_indicate, phy))

phy_btn.place(x=15,y=400)

phy_indicate = tk.Label(frame, text='', bg='#1E3694')
phy_indicate.place(x=3,y=400, width=5, height=85)

def phy():
    choose_file = tk.Label(main_frame,text="Export a fasta file",bg='#5B72B4', fg = 'Black', font=('Bold', 15))
    choose_file.place(x=10,y=30)

    writeseq = tk.Label(main_frame,text="or write your sequences here",bg='#5B72B4', fg = 'Black', font=('Bold', 10))
    writeseq.place(x=10,y=75)

    align = tk.Button(main_frame,text='Save input',command=lambda:seq_tree(), font=('Bold', 10),
                        fg='White', bd=0, bg='#91A1CD')
    align.place(x=250,y=70)

    file_btn = tk.Button(main_frame, text='Browse', font=('Bold', 10), width=20,
                        fg='White', bd=0, bg='#91A1CD', command=lambda:browseFiles(3))
    file_btn.place(x=200,y=35)

    UPGMA_btn = tk.Button(main_frame, text='UPGMA', font=('Bold', 20),
                        fg='White', bd=0, bg='#1E3694', command=lambda:do_phy())
    UPGMA_btn.place(x=1000,y=200)

    Newick_btn = tk.Button(main_frame, text='Newick', font=('Bold', 20),
                        fg='White', bd=0, bg='#1E3694', command=lambda:do_dnd())
    Newick_btn.place(x=1000,y=400)

    global outputtree
    outputtree = tk.Text(main_frame, height = 25,
              width = 85,
              bg = "#91A1CD")
        
    outputtree.place(x=30,y=100)








# Statistics
stat_btn = tk.Button(frame, text='Statistics', font=('Bold', 20),
                        fg='White', bd=0, bg='#91A1CD', command=lambda:indicate(stat_indicate,stat))

stat_btn.place(x=40,y=530)

stat_indicate = tk.Label(frame, text='', bg='#1E3694')
stat_indicate.place(x=3,y=530, width=5, height=55)

def stat():
    
    choose_file = tk.Label(main_frame,text="Export a fasta file",bg='#5B72B4', fg = 'Black', font=('Bold', 20))
    choose_file.place(x=450,y=30)

    file_btn = tk.Button(main_frame, text='Browse', font=('Bold', 12), width=20,
                        fg='White', bd=0, bg='#91A1CD', command=lambda:browseFiles(2))
    file_btn.place(x=700,y=35)

    writeseq = tk.Label(main_frame,text="or write your sequence here",bg='#5B72B4', fg = 'Black', font=('Bold', 10))
    writeseq.place(x=450,y=75)

    process = tk.Button(main_frame,text='Done!',command=lambda:seq_stat(), font=('Bold', 10),
                        fg='White', bd=0, bg='#91A1CD')
    process.place(x=700,y=75)

    global seqoutput
    seqoutput = tk.Text(main_frame, height = 10,
              width = 50,
              bg = "#91A1CD")
        
    seqoutput.place(x=400,y=100)


    lengthlbl = tk.Label(main_frame,text="~Sequence Length:",bg='#5B72B4', fg = 'Black', font=('Bold', 15))
    lengthlbl.place(x=100,y=400)

    baseslbl = tk.Label(main_frame,text="~Base composition:",bg='#5B72B4', fg = 'Black', font=('Bold', 15))
    baseslbl.place(x=100,y=500)

    gclbl = tk.Label(main_frame,text="~GC Content:",bg='#5B72B4', fg = 'Black', font=('Bold', 15))
    gclbl.place(x=100,y=600)
 
    global length
    global bases
    global gc

    bases = tk.Text(main_frame, height = 2,
              width = 30,
              bg = "#91A1CD")
        
    bases.place(x=350,y=500)

    length = tk.Text(main_frame, height = 2,
              width = 10,
              bg = "#91A1CD")
        
    length.place(x=350,y=400)

    gc = tk.Text(main_frame, height = 2,
              width = 20,
              bg = "#91A1CD")
        
    gc.place(x=350,y=600)




















frame.pack(side=tk.LEFT)
frame.pack_propagate(False)
frame.configure(width=200, height=700)

main_frame = tk.Frame(root, bg='#5B72B4' ,highlightbackground='white',
                        highlightthickness=3)

main_frame.pack(side=tk.LEFT)
main_frame.pack_propagate(False)
main_frame.configure(width=1500, height=700)


root.mainloop()