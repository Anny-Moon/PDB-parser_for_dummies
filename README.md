# PDB Parser for Dummies

Here is a very simple reader for [Protein Data Bank (PDB)](http://www.rcsb.org) files.
(At this stage of the project) Reader can extract xyz-coordinates of C-alpha (CA)
atoms from the first (or the only one) model and writes the result into
dat-file in format:

```
<x_firts_atom>    <y_first_atom>    <z_first_atom>
<x_second_atom>   <y_second_atom>   <z_second_atom>

 ...               ...               ...
 
<x_last_atom>     <y_last_atom>     <z_last_atom>
```
If there are any missing atoms in the model, the program will tell about it
and show maps for them both in percents and actual size.

You will be able to rewrite dat-file with only one segment in order to not have any
missings.

If there are any CA atoms with the same number, the program automatically
takes only the first and ignore all the rest with the same number.

##Example: 5dn7
```
$ ./pdb-reader 5dn7
```
The output is:
```
Protein: 5dn7
Missing atoms from 361 to 365 (5 atoms) 
Missing atoms from 558 to 562 (5 atoms) 
Number of C-alpha atoms in the model: 250
But there is data only for 240 of them.

*****************************************
* Maps of missing atoms                 *
* atom: .         missing atom: 0       *
*****************************************
Percentage: string length = 100 chars.

......00.............................................................................00.............

Actual: string length = number of atoms in model.

................00000................................................................................................................................................................................................00000................................

*****************************************
If you want to rewrite dat-file with only one segment call
./pdb-reader (without arguments) for instructions.
```

The program created file `results/xyz_5dn7.dat` with 240 lines.

Let's concider the maps. They represent the positions and amounts of missing atoms, where `0` stands for missing atom and `.` for atom from the pdb-file.
The first map depicture the percentege picture (so not very precise becouse of rounding). The second one is
an actual picture: each char is one char.

Now I want to rewrite the dat-file with one segment in oder to not have missing data there. There are 3 segments
there. I will take the largest one, between two missing parts. I can do this by calling the program again but
now give the number of the first atom (or any missing atom before the first) in segment as the **second**
argument.

Like this:
```
$ ./pdb-reader 5dn7 365
```
or the same:
```
$ ./pdb-reader 5dn7 366
```
## For Whom It Can Be Useful

 

Anna Sinelnikova
Uppsala, 2017
