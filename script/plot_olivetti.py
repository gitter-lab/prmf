import sys, argparse
from sklearn.datasets import fetch_olivetti_faces
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="""
  Compose one image containing all faces in the Olivetti faces dataset
  corresponding to a single person, identified by <--person-id>.
""")
  parser.add_argument("--person-id", default=0, type=int)
  parser.add_argument("outfile")
  args = parser.parse_args()

  dataset = fetch_olivetti_faces()
  image_vecs = dataset.data
  person_ids = dataset.target

  inds = np.where(person_ids == args.person_id)
  if(len(inds) > 0):
    select_vecs = image_vecs[inds]
    gallery_col = np.ceil(np.sqrt(len(select_vecs)))
    gallery_row = gallery_col
    plt.figure(figsize=(2. * gallery_col, 2.26 * gallery_row))
    for i, image_vec in enumerate(select_vecs):
      plt.subplot(gallery_row, gallery_col, i+1)
      plt.imshow(image_vec.reshape((64,64)), cmap=plt.cm.gray)
    plt.savefig(args.outfile)
  else:
    sys.stderr.write("invalid person_id: {}\n".format(args.person_id))
    sys.exit(21)
