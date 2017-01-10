# PDB Parser for Dummies

Here is a very simple reader for [Protein Data Bank (PDB)](http://www.rcsb.org) files.
(At this stage of the project) Reader can extract xyz-coordinates of C-alpha (CA)
atoms from the first (or the only one) model and writes the result into
dat-file in format:

<x_firts_atom>    <y_first_atom>    <z_first_atom>
 ...               ...               ...
<x_last_atom>     <y_last_atom>     <z_last_atom>

If there are any missing atoms in the model, the program will tell about it
and show maps for them both in percents and actual size.

You will be able to rewrite dat-file with only one segment in order to not have any
missings.

If there are any CA atoms with the same number, the program automatically
takes only the first and ignore all the rest with the same number.

For exmple, for protein 5dn7 the percentage map is following:

......00.............................................................................00.............

The length of the line is 100 chars.
## For Whom It Can Be Useful

 

Anna Sinelnikova
Uppsala, 2017