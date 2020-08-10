
# SDG tutorial

[TOC]

# Note: Scope of this tutorial

This tutorial is meant to reflect current practices to use SDG. It will track "a current version of SDG", and as such won't keep compatibility with the versions released in publications.

# 1. The SDG Workspace

The SDG 

```python

```



# 2. Working with a SequenceDistanceGraph

## Loading and exporting graphs

## SDG graph representation

## Navigating the graph: NodeViews and LinkViews

# 3. K-mer frequencies

# 4. Datastores and mappers

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