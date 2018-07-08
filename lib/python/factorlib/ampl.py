def ampl_write_sparse_arrs(arrs, out_fh, arr_name="L", set_name="L_set", default=0):
  """
  Write to an AMPL data file Ls in the following form:

  param L default 0 := 
    1 1 1 1
    1 1 2 -1
    ...
    ;

  where first dimension is the array index, second dimension is row, third dimension is column

  Parameters
  ----------
  arrs : list of scipy Compressed Sparse Row matrix

  out_fh : io-like

  arr_name : str
  
  default : float
    default value for not array entries not mentioned in sparse arrays
  """
  arrs = list(arrs)
  indent = "  "
  set_vs = " ".join(map(str, map(lambda x: x+1, range(len(arrs)))))
  out_fh.write("set {} := {};\n".format(set_name, set_vs))
  out_fh.write("param {} default {} :=\n".format(arr_name, default))
  for i, arr in enumerate(arrs):
    i_arr, j_arr = arr.nonzero()
    for j, k in zip(i_arr, j_arr):
      v = arr[j,k]
      out_fh.write(indent + " ".join(map(str, [i+1,j+1,k+1,v])) + "\n")
  out_fh.write(";")

def ampl_write_arr_shape(out_fh, n_row, n_col, row_name="I", col_name="J"):
  row_str = " ".join(map(str, range(1, n_row+1)))
  col_str = " ".join(map(str, range(1, n_col+1)))
  out_fh.write("set {} := {};\n".format(row_name, row_str))
  out_fh.write("set {} := {};\n".format(col_name, col_str))

def ampl_write_sparse_arr(arr, out_fh, n_nodes, index_map=None, arr_name="X", row_name="I", col_name="J", default=0):
  """
  Single array version of ampl_write_sparse_arrs (that is, the matrix verison, rather than 
  the 3D tensor version). However, the use case is more to embed a dense block matrix inside
  a larger sparse matrix (which is only sparse because the entries outside the block are 0).
  This situation arises when we have measurements on a limited set of genes but wish to
  define a model on

  Parameters
  ----------
  arr

  out_fh

  index_map : dict or None
    if dict then mapping of array index to array index where the domain index is (origin-0) for <arr>
    and the range index is (origin-0) for a larger vector space
  """
  i_arr, j_arr = arr.nonzero()
  n_rows, n_cols = arr.shape
  ampl_write_arr_shape(out_fh, n_nodes, n_cols)
  out_fh.write("param {} default {} :=\n".format(arr_name, default))
  for i, j in zip(i_arr, j_arr):
    v = arr[i,j]
    i_map = j_map = None
    if index_map is not None:
      i_map = index_map[i]
      j_map = index_map[j]
    else:
      i_map = i
      j_map = j
    out_fh.write(" ".join(map(str, [i_map+1, j_map+1, v])) + "\n")

def write_ampl_laplacians(Ls, out_fh, arr_name="L"):
  """
  Write to an AMPL data file Ls in the following form:

  param L := 
    [1, *, *]: 1 2 3 :=
      1 1 -1 0
      2 -1 2 -1
      3 0 -1 1
    [2, *, *]: 1 2 3 :=
      1 1 -1 0
      2 -1 1 0
      3 0 0 0
    ;
  """
  indent = "  "
  out_fh.write("param {} :=\n".format(arr_name))
  for i in range(len(Ls)):
    L = Ls[i]
    n_row, n_col = L.shape
    col_str = " ".join(map(str, range(1,n_col+1)))
    out_fh.write(indent + "[{}, *, *]: {} :=\n".format(i+1, col_str))
    for j in range(n_row):
      out_fh.write(indent*2 + " ".join(map(str, [j+1] + list(L[j,:]))) + "\n")
  out_fh.write(indent + ";")

def write_ampl_data(data, out_fh, data_name="X", row_name="I", col_name="J"):
  """
  Write an ampl data file similar to the following example:

  set I := 1 2 3;
  set J := 1 2 3;
  param X: 1 2 3 :=
    1 200 100 150
    2 10 15 200
    3 100 50 0 ;
  """
  n_row, n_col = data.shape

  row_str = " ".join(map(str, range(1, n_row+1)))
  col_str = " ".join(map(str, range(1, n_col+1)))
  #"1 .. {}".format(n_row)
  #"1 .. {}".format(n_col)

  out_fh.write("set {} := {};\n".format(row_name, row_str))
  out_fh.write("set {} := {};\n".format(col_name, col_str))
  out_fh.write("param {}: {} :=\n".format(data_name, col_str))

  for i in range(n_row):
    row_num = i+1
    out_fh.write("  " + str(row_num) + " " + " ".join(map(str, data[i,:])) + "\n")
  out_fh.write(";\n")

def write_ampl_params(n_components, out_fh, comp_name="K_card"):
  """
  Write an ampl data file in the following form:

  param K_card := 2;
  """
  out_fh.write("param {} := {};".format(comp_name, n_components))
