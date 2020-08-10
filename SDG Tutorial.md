# SDG tutorial

SDG tutorial
[TOC]

# Note: Scope of this tutorial

This tutorial is meant to reflect current practices to use SDG. It will track "a current version of SDG", and as such won't keep compatibility with the versions released in publications.

# 1. The SDG Workspace

The SDG workspace...

```python
>>> import SDGpython as SDG

>>> ws = SDG.WorkSpace()
```

now you can acess a workspace using the ws variable.

To check the status of your workspace you can use the `ws.ls()` method

``` python
>>> print(ws.ls())
SDG Workspace
  SequenceDistanceGraph SDG: 0 nodes, 0 links
  Distance Graphs: 0
  Paired Read Datastores: 0
  Linked Read Datastores: 0
  Long Read Datastores: 0
  Kmer Count Datastores: 0
```

at this point the workspace is empty.


# 2. Working with a SequenceDistanceGraph 

Explanation of what a SDG is

## SDG graph representation

Description of a Sequence distance graph how its implemented on SDG

## Loading and exporting graphs

SDG accetps gfa1 graphs as long as bcalm graphs and plain fasta files (fasta files are loaded as disconected graphs). To add a graph to an existing workspace

```python
## Loading an ecoli graph created using sdg-dbg
>>> ws.sdg.load_from_gfa('./ecoli_assm_DBG.gfa')
2020-08-10 16:06:50: Graph fasta filesname: ./ecoli_assm_DBG.fasta
2020-08-10 16:06:50: Loading sequences from ./ecoli_assm_DBG.fasta
2020-08-10 16:06:50: 811 nodes loaded (0 canonised).
2020-08-10 16:06:50: 811 nodes after connecting with 1091 links.
```

Now the graph is contained within the workspace, if we run now the report we'll observer the available graph in the report

```python
>>> print(ws.ls())
SDG Workspace
  SequenceDistanceGraph SDG: 811 nodes, 1092 links
  Distance Graphs: 0
  Paired Read Datastores: 0
  Linked Read Datastores: 0
  Long Read Datastores: 0
  Kmer Count Datastores: 0
  
```
The main sequence distance graph of the workspace can be referenced by acesssing `ws.sdg`.
This is the main graph of the workspace and stores the genome sequence and the links that join the sequences.
Besides the `SequenceDistanceGraph` (SDG) in sdg you can create `DistanceGraph` (DG) the main difference between an SDG and a DG is that the DG only stores links between nodes referencing the sequence directly from the main SDG of the workspace.
DG are used to express a different linkage of the original graph, usually this is done during the assembly process to store the result of some graph processing/ordering.

To write a graph to disk 

```python
>>> ## If you want to output the main graph
>>> ws.sdg.write_to_gfa1('./output_graph_name.gfa')
```

as SDG are a special version of DG (DG+sequence) the SDG implementation within the SDG library inherits all methods from DG. that means that most of the SDG methid can be used with DG and vice versa with some exceptions.

## Navigating the graph: NodeViews and LinkViews

To explore the graph SDG has a structure called Views, si far there are 3 implemented types of views. `NodeView` (centered on nodes), `LinkView` (centered on links) and `TangleView` (centered on complex regions of the graph).

The most usefull is the `NodeView`, nodeviews "sit" on top of nodes to acess a number of possible attributes of the particular node you are "viewing".
For example, to acess node #10

```python
>>> nv = ws.sdg.get_nodeview(10)
>>> print(nv)
<NodeView: node 10 in graph SDG>
```

now using `nv` we can acess a number of node characteristics, for example

```python
## Node id
>>> nv.node_id()
10

## Size
>>> nv.size()
103

## Sequence
>>> nv.sequence()
'GCCTGATACAGATTAAATCAGAACGCAGAAGCGGTCTGATAAAACAGAATTTGCCTGGCGGCCGTAGCGCGGTGGTCCCACCTGACCCCATGCCGAACTCAGA'

## Reverse complement the node and then get the sequence
>>> nv.rc().sequence()
'TCTGAGTTCGGCATGGGGTCAGGTGGGACCACCGCGCTACGGCCGCCAGGCAAATTCTGTTTTATCAGACCGCTTCTGCGTTCTGATTTAATCTGTATCAGGC'

## Next node in the graph is node 27 overlaping 62 bases (k-1)
>>> nv.next()
[<LinkView to 27 at -62 bp>]

## Prev node in the graph is node 323 overlapping 62 bases (k-1)
>>> nv.prev()
[<LinkView to -323 at -62 bp>]

```
This is just a sample of all available methods for a nodeview. The available methis depend on the available workspace information and include kmer coverage, mapped reads to the node, tags on the node if there is linked read information, etc.
 More method are going to appear as the different sections are explained.
 
 for a graph all nodeviews can be retrieved using the `ws.sdg.get_all_nodeviews()` 
 
```python
node_sizes = []
for nv in ws.sdg.get_all_nodeviews():
    node_sizes.append(nv.size())
node_sizes.sort()
```
this returns a vector with all nodeviews for the graph, each element of the vector is a nodeview for a corresponding node in the graph.

 
 
# 3. K-mer frequencies

Workspaces can hold multiple kmer collections inside, kmers are stored in `kmer_counter`s containers inside the workspace.

```python
kmer_counter = ws.add_kmer_counter('pe_k31', 31)
```

You can check the kmer counter by doing `ws.ls()`

```python
>>> print(ws.ls())
SDG Workspace
  SequenceDistanceGraph SDG: 811 nodes, 1092 links
  Distance Graphs: 0
  Paired Read Datastores: 0
  Linked Read Datastores: 0
  Long Read Datastores: 0
  Kmer Count Datastores: 0
    KmerCounter pe_k31: index with 4555914 31-mers
      Counts: 1
        sdg: 4595583 total 31-mers
```
When a kmer counter is created a first count is added, this initial count stores the graph kmers at the specified k value (31 in this case). This count reflects how many times each 31-mer of the graph is present on that same graph. This feature will come usefull when we analyze repetitions, heterozygosity or work with kci values in any way.

Then extra counts can be added to represent different datasets

```python
>>> kmer_counter.add_count('pe', ['./pe/ecoli_pe_R1.fastq.gz', './pe/ecoli_pe_R2.fastq.gz'])
2020-08-10 19:44:05: Populating lookup map
2020-08-10 19:44:07: Map populated with 4555914 entries
2020-08-10 19:44:07: Counting from file: ./pe/ecoli_pe_R1.fastq.gz
2020-08-10 19:44:09: 99448 reads processed 23999992 / 26781064 kmers found
2020-08-10 19:44:12: 198152 reads processed 47999992 / 53364137 kmers found
2020-08-10 19:44:14: 296331 reads processed 71999992 / 79799347 kmers found
2020-08-10 19:44:16: 394409 reads processed 95999992 / 106215960 kmers found
2020-08-10 19:44:17: 500000 reads processed 121676244 / 134649543 kmers found
2020-08-10 19:44:17: Counting from file: ./pe/ecoli_pe_R2.fastq.gz
2020-08-10 19:44:20: 590481 reads processed 139676236 / 159022475 kmers found
2020-08-10 19:44:22: 699698 reads processed 161676236 / 188438406 kmers found
2020-08-10 19:44:23: 798314 reads processed 181676236 / 214998725 kmers found
2020-08-10 19:44:26: 896999 reads processed 201676236 / 241583986 kmers found
2020-08-10 19:44:27: 1000000 reads processed 222481360 / 269324800 kmers found
2020-08-10 19:44:27: Done

```

The kmers will be stored inside the kmer counter in a kmer collection named 'pe'.
Once the kmers are counted and stored in the ws we can get some basic stas using the `ws.ls()` workspace function.

```python
>>> print(ws.ls())
SDG Workspace
  SequenceDistanceGraph SDG: 811 nodes, 1092 links
  Distance Graphs: 0
  Paired Read Datastores: 0
  Linked Read Datastores: 0
  Long Read Datastores: 0
  Kmer Count Datastores: 0
    KmerCounter pe_k31: index with 4555914 31-mers
      Counts: 2
        sdg: 4595583 total 31-mers
        pe: 222481360 total 31-mers
```

Now the new kmer count is shown in the pe_k31 collection.
Now add a second count from the pb reads k31 also.

```python
>>> kmer_counter.add_count('pb', ['./pb/ecoli_pb_all.fastq',])
2020-08-10 19:57:11: Populating lookup map
2020-08-10 19:57:12: Map populated with 4555914 entries
2020-08-10 19:57:12: Counting from file: ./pb/ecoli_pb_all.fastq
2020-08-10 19:57:33: 91394 reads processed 14432912 / 479861967 kmers found
2020-08-10 19:57:33: Done

>>> print(ws.ls())
SDG Workspace
  SequenceDistanceGraph SDG: 811 nodes, 1092 links
  Distance Graphs: 0
  Paired Read Datastores: 0
  Linked Read Datastores: 0
  Long Read Datastores: 0
  Kmer Count Datastores: 0
    KmerCounter pe_k31: index with 4555914 31-mers
      Counts: 3
        sdg: 4595583 total 31-mers
        pe: 222481360 total 31-mers
        pb: 14432912 total 31-mers

```

Now the pacbio reads 31-mers are also stored in the same counter along the pe count.

```python
>>> import matplotlib.pylab as plt
>>> plt.plot(nv.kmer_coverage('pe_k31', 'pe'))
>>> plt.ylabel("# kmers")
>>> plt.xlabel("bp")
>>> plt.show()

```

![pe 31-mers coverage](https://i.imgur.com/9nqrk7r.png)


Other kmer collections at other kmer values can also be added to the workspace, for example, pacbio reads kmer at k=17.

```python
>>> kmer_counter_pb = ws.add_kmer_counter('pb_k17', 17)
>>> kmer_counter_pb.add_count('pb', ['./pb/ecoli_pb_all.fastq',])
2020-08-10 20:46:05: Populating lookup map
2020-08-10 20:46:07: Map populated with 4531666 entries
2020-08-10 20:46:07: Counting from file: ./pb/ecoli_pb_all.fastq
2020-08-10 20:46:35: 91394 reads processed 60543658 / 481141483 kmers found
2020-08-10 20:46:35: Done

>>> print(ws.ls())
SDG Workspace
  SequenceDistanceGraph SDG: 811 nodes, 1092 links
  Distance Graphs: 0
  Paired Read Datastores: 0
  Linked Read Datastores: 0
  Long Read Datastores: 0
  Kmer Count Datastores: 0
    KmerCounter pe_k31: index with 4555914 31-mers
      Counts: 3
        sdg: 4595583 total 31-mers
        pe: 222481360 total 31-mers
        pb: 14432912 total 31-mers
    KmerCounter pb_k17: index with 4531666 17-mers
      Counts: 2
        sdg: 4606937 total 17-mers
        pb: 60543658 total 17-mers
```

Now there is a new kmer counter collection at k17. with it's corresponding sdg counter and the pacbio kmers collection.


# 4. Datastores and mappers

SDG supports 3 types of generic data types. Paired reads (pe, mp), Linked reads (10x type) and long reads (pacbio, nanopore). 

# 5. Disecting a graph: KCI and Tangles and Tips

# 6. Linkage analyses: using DistanceGraphs

# 7. Editing the graph: the GraphEditor

# 8. Common graph simplification algorithms: the GraphContigger

# 9. Untangling and reconnecting: strider.

```
t=c=0
while ws.sdg.get_all_nodeviews():
    if len(nv.prev())==1 and len(nv.next())==1 and len(nv.parallels())==1:
        t+=nv.size()
        c+=1
print(t,c)
```



<!--stackedit_data:
eyJoaXN0b3J5IjpbMjM3MTQ3OTczXX0=
-->