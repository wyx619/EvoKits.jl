# GoM reference

`reference.tsv` contains three tiny fixed configurations evaluated with local
maptpx 1.9-7 under R 4.6.2: the scalar-shape NEF log objective and one native
`tpxEM` update, including an all-zero species column and a three-region model.
Matrices are packed in column-major order. Values use 17 significant digits.
The unnormalized objective uses `shape*sum(log(theta)) + sum(log(omega))/k`;
`shape` is not reduced by one. Ordinary tests load this table without R.

Source: the user-provided maptpx 1.9-7 source distribution, GPL-3. The local
generator and larger fit comparisons are kept under ignored `validation/`.
