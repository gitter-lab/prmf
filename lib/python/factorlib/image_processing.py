OLIVETTI_SHAPE = (64, 64)

def plot_olivetti_components(plt, images, title='Components'):
  """
  Parameters
  ----------
  images : list of np.array
  """
  gallery_row = 2
  gallery_col = 3
  image_shape = OLIVETTI_SHAPE
  plt.figure(figsize=(2. * gallery_col, 2.26 * gallery_row))
  plt.suptitle(title, size=16)
  for i, comp in enumerate(images):
    plt.subplot(gallery_row, gallery_col, i + 1)
    vmax = max(comp.max(), -comp.min())
    plt.imshow(comp.reshape(image_shape), cmap=plt.cm.gray,
           interpolation='nearest',
           vmin=-vmax, vmax=vmax)
    plt.xticks(())
    plt.yticks(())
  plt.subplots_adjust(0.01, 0.05, 0.99, 0.93, 0.04, 0.)

def pixel_similarity(n_row, n_col, i, j):
  x_i = i // n_row
  y_i = i % n_row
  x_j = j // n_row
  y_j = j % n_row
  dist = math.sqrt((x_i - x_j) ** 2 + (y_i - y_j) ** 2)
  return dist

def pixel_to_index(n_col, i, j):
  # Convert from pixel coordinates to vector coordinates
  return i * n_col + j

def pixel_prior(mat):
  """
  In images, we expect the pixel intensity of nearby pixels to be only slightly different 
  from the neighboring pixels
  """
  G_pixel = nx.Graph()
  n_row, n_col = mat.shape
  for i,j in it.product(range(n_row), range(n_col)):
    pixel_v = mat[i,j]
    neighbors = get_neighbors(n_row, n_col, i,j)

    for neighbor in neighbors:
      k_1 = pixel_to_index(n_col, i,j)
      k_2 = pixel_to_index(n_col, *neighbor)

      G_pixel.add_edge(k_1, k_2, {'weight': 1})
  return G_pixel

def pixel_nn_prior(mat, percentile=50):
  """
  Identify the neighboring pixel with the most similar intensity
  """
  G_nn = nx.Graph()
  n_row, n_col = mat.shape
  vec = mat.flatten()
  percentile_v = np.percentile(vec, percentile)
  stdev = np.std(vec)
  for i,j in it.product(range(n_row), range(n_col)):
    pixel_v = mat[i,j]
    if(pixel_v < percentile_v):
      continue

    neighbors = get_neighbors(n_row, n_col, i,j)

    best = None
    best_sim = 0
    for neighbor in neighbors:
      neigh_v = mat[neighbor[0], neighbor[1]]
      sim = rbf_similarity(stdev, pixel_v, neigh_v)
      if neigh_v >  percentile_v and sim > best_sim:
        best = neighbor
        best_sim = sim

    if best is not None:
      k_1 = pixel_to_index(n_col, i,j)
      k_2 = pixel_to_index(n_col, *best)
      G_nn.add_edge(k_1, k_2, {'weight': best_sim})

  return G_nn

def vec_to_graph_comb(vec, percentile=50, n_neighbors=1, n_row=64, n_col=64):
  """
  Construct network that is a combination of nearest_neighbors in the image space
  and pixel intensity similarity among those neighbors
  """
  mat = vec.reshape((n_row, n_col))
  G_pixel = pixel_prior(mat)
  G_nn = pixel_nn_prior(mat, percentile)
  G = combine_graphs([G_pixel, G_nn], [0.5, 0.5])
  return G

def vec_to_graph_prim(vec, percentile=50, n_neighbors=1, n_row=64, n_col=64):
  """
  Primitive version which only examines the neighboring pixels; always 1-NN
  """
  G = nx.Graph()

  def pixel_to_index(i, j):
    # Convert from pixel coordinates to vector coordinates
    return i * n_col + j

  stdev = np.std(vec)
  mat = vec.reshape((n_row, n_col))
  percentile_v = np.percentile(vec, percentile)

  for i,j in it.product(range(n_row), range(n_col)):
    neighbors = []
    pixel_v = mat[i,j]
    if(pixel_v < percentile_v):
      continue

    if i != 0:
      neighbors.append((i-1, j))
    if i != n_row - 1:
      neighbors.append((i+1, j))
    if j != 0:
      neighbors.append((i, j-1))
    if j != n_col - 1:
      neighbors.append((i, j+1))

    best = None
    best_dist = 0
    for neighbor in neighbors:
      neigh_v = mat[neighbor[0], neighbor[1]]
      dist = rbf_similarity(stdev, pixel_v, neigh_v)
      if neigh_v >  percentile_v and dist > best_dist:
        best = neighbor
        best_dist = dist

    if best is not None:
      k_1 = pixel_to_index(i,j)
      k_2 = pixel_to_index(*best)
      G.add_edge(k_1, k_2, {'weight': best_dist})

  return G

def get_neighbors(n_row, n_col, i,j):
  """
  Get neighboring pixel coordinates
  """
  neighbors = []
  if i != 0:
    neighbors.append((i-1, j))
  if i != n_row - 1:
    neighbors.append((i+1, j))
  if j != 0:
    neighbors.append((i, j-1))
  if j != n_col - 1:
    neighbors.append((i, j+1))
  return neighbors

# TODO n_row and n_col 
def vec_to_graph(vec, percentile=50, n_neighbors=1, n_row=64, n_col=64):
  vec = vec.flatten()
  stdev = np.std(vec)
  percentile_v = np.percentile(vec, percentile)

  def rbf_sim(x,y):
    return rbf_similarity(stdev, x, y)

  def sim(vec, i, j, alpha=0.5):
    return alpha * rbf_sim(vec[i], vec[j]) + (1-alpha) * pixel_similarity(n_row, n_col, i, j)

  G = nx.Graph()
  # TODO slow implementation
  for i, j in it.combinations(range(vec.shape[0]),2):
    # <percentile> percent edge weight threshold
    if(vec[i] > percentile_v and vec[j] > percentile_v):
      G.add_edge(i,j, {'weight': sim(vec, i, j, alpha=0.1)})
  G_prime = transform_nearest_neighbors(G, n_neighbors)
  return G_prime
