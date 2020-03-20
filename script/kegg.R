#!/usr/bin/env Rscript
load_success = TRUE
load_success = load_success && library(KEGGREST, quietly=TRUE, logical.return=TRUE)
load_success = load_success && library(KEGGgraph, quietly=TRUE, logical.return=TRUE)
load_success = load_success && library(graph, quietly=TRUE, logical.return=TRUE)
load_success = load_success && library(igraph, quietly=TRUE, logical.return=TRUE)
load_success = load_success && library(argparse, quietly=TRUE, logical.return=TRUE)
if(!load_success) {
  write("Unable to load dependencies", stderr())
  quit(status=2)
}

HSA_PREFIX = "hsa:"
# [1] "NCBI-ProteinID: NP_002761"     "NCBI-GeneID: 5645"
# [3] "OMIM: 601564"                  "HGNC: 9483"
# [5] "HPRD: 03335"                   "Ensembl: ENSG00000275896"
# [7] "Vega: OTTHUMG00000158904"      "Pharos: P07478(Tchem)"
# [9] "UniProt: P07478 Q5NV56 Q6PK75" 
ENSEMBL = "Ensembl"
HGNC = "HGNC"
REFSEQ = "NCBI-ProteinID"
PATHWAY_DELIM = ":"
LINK_DELIM = ":"

#' Return a graph object associated with a kegg pathway identifier.
#' Parse existing KGML file if it exists. Otherwise download, save, then parse.
#'
#' @param pathway_id e.g. "hsa:05164"
#' @return graph object (see help(graphNEL))
get_pathway <- function(pathway_id, outdir=NULL) {
  # remove ":" if present
  pathway_id = sub(':', '', pathway_id)
  pathway_file = NULL
  if (is.null(outdir)) {
    pathway_file = tempfile()
  } else {
    pathway_bn = paste(pathway_id, '.xml', sep='')
    pathway_file = file.path(outdir, pathway_bn)
  }
  kegg.pathway = NULL
  if (!file.exists(pathway_file)) {
    pathway_org = substr(pathway_id, 1, 3)
    pathway_num = substr(pathway_id, 4, nchar(pathway_id))
    kgml = retrieveKGML(pathwayid=pathway_num, organism=pathway_org, destfile=pathway_file)
  }
  kegg.pathway = parseKGML(pathway_file)
  graphnel = KEGGpathway2Graph(kegg.pathway, expandGenes = TRUE)
  graph = graph_from_graphnel(graphnel)
  return(graph)
}

#' Translate KEGG node identifiers using id_map. id_map is a one-to-many mapping. Create extra nodes
#' as required.
#'
#' TODO figure out name of 'weight' edge attribute from G rather than hard code
rename_nodes <- function(G, id_map) {
  nodes = V(G)
  # first, create new nodes for those mapped to >= 2 ids where rename mapping
  # is done implicitly through node creation
  for(node_ind in nodes) {
    node = nodes[node_ind]
    old_name = node$name
    new_names = id_map[[old_name]]
    if (length(new_names) >= 2) {
      # new nodes
      for (i in 2:length(new_names)) {
        new_name = new_names[i]
        G = add_vertices(G, 1, name=new_name)
      }

      # "from" edges
      edge_sequence = E(G)[from(old_name)]
      ends_rv = ends(G, edge_sequence)
      if (length(edge_sequence) > 0) {
        for (i in 2:length(new_names)) {
          new_name = new_names[i]
          for (j in 1:length(edge_sequence)) {
            # j indexes which edge in edge sequence, 2 indexes the target of the edge
            target = ends_rv[j,2]
            G = add_edges(G, c(new_name, target), weight=1)
          }
        }
      }

      # "to" edges
      edge_sequence = E(G)[to(old_name)]
      ends_rv = ends(G, edge_sequence)
      if (length(edge_sequence) > 0) {
        for (i in 2:length(new_names)) {
          new_name = new_names[i]
          for (j in 1:length(edge_sequence)) {
            # j indexes which edge in edge sequence, 1 indexes the source of the edge
            source_node = ends_rv[j,1]
            G = add_edges(G, c(source_node, new_name), weight=1)
          }
        }
      }
    }
  }

  # then, rename nodes
  for(node_ind in nodes) {
    node = nodes[node_ind]
    old_name = node$name
    new_names = id_map[[old_name]]
    if (length(new_names) >= 1) {
      G = set_vertex_attr(G, 'name', index=old_name, value=new_names[1])
    }
  }
  return(G)
}

#' Write graph object to a file
#'
#' @param 
#' @todo try igraph package because graph does not have its own serialization?  #' @todo what about singleton nodes?
#' @todo this will write extra columns when an hsa maps to multiple ensemble identifiers (which occurs in hsa:05164)
write_abc <- function(graph, outfile, id_map = NA) {
  fh <- file(outfile, open="w")

  write_edge <- function(src, target) {
    line = sprintf("%s\t%s", src, target)
    writeLines(line, fh)
  }
  write_edges <- function(src, targets) {
    if (length(targets) > 0) {
      mapply(write_edge, src, targets)
    }
  }

  edges_rv = edges(graph)
  srcs = names(edges_rv) # src node labels
  targets = edges_rv # adjacency list like representation of edges

  # relabel nodes if id_map was provided
  if(is.na(id_map)) {
    # then do not relabel nodes
  } else {
    # then we use the id_map to relabel nodes
    srcs_new = mapply(function(src) {
      # TODO if src is missing and value is NULL, how does mapply treat this
      return(id_map[[src]])
    }, srcs)
    targets_new = mapply(function(target_list) {
      target_list_new = mapply(function(target) {
        return(id_map[[target]])
      }, target_list)
      return(target_list_new)
    }, targets)

    srcs = srcs_new
    targets = targets_new
  }

  mapply(write_edges, src=srcs, targets=targets)
  close(fh)
}

# TODO number of genes off by 1 when compared to get_pathway
genes_in_pathway_http <- function(pathway_id) {
  pathway_id = gsub(PATHWAY_DELIM, '', pathway_id)
  pathway_data = keggGet(c(pathway_id))

  # odd entries are homo sapien identifiers
  genes = pathway_data[[1]]$GENE
  gene_ids = array(dim=floor((length(genes)-1)/2))
  for(i in 1:length(genes)) {
    if(i %% 2 == 1) {
      # then odd
      gene_id = paste(HSA_PREFIX, genes[i], sep="")
      gene_ids[floor(i/2)] = gene_id
    }
  }

  return(gene_ids)
}

genes_in_pathway <- function(graph) {
  return(vertex_attr(graph, 'name', V(graph)))
}

#' Link KEGG gene identifiers e.g. "hsa:5291" to the identifier named in <type>, default ENSEMBL
#'
#' @return ids_map an "environment" used as a hash mapping KEGG gene ids to their ENSEMBL counterparts
#'   (one-to-many) or an empty array if there is no linked identifier
link_genes_hash <- function(kegg_gene_ids, type=ENSEMBL) {
  linked_genes = link_genes_array(kegg_gene_ids, type)
  ids_map = new.env(parent=emptyenv())
  mapply(function(id, value) {
    ids_map[[id]] = value
  }, id=kegg_gene_ids, value=linked_genes)
  return(ids_map)
}

#' Link KEGG gene identifiers e.g. "hsa:5291" to the identifier named in <type>, default ENSEMBL
#'
#' @return link_ids an array of the same length as kegg_gene_ids with the mapped identifiers 
#'   (note one-to-many mapping) or an empty array if there is no linked identifier
link_genes_array <- function(kegg_gene_ids, type=ENSEMBL) {
  link_ids = array(dim = length(kegg_gene_ids))

  # return character vector of length 2: <id_type>, <id_value>
  # or an empty vector if link cannot be parsed
  parse_link = function(dblink) {
    rv = c()
    split_rv = strsplit(dblink, ':', fixed = TRUE)
    if(split_rv == dblink) {
      # then split failed, return NULL
    } else {
      id_type = trimws(split_rv[[1]][1])
      id_value_str = trimws(split_rv[[1]][2])
      id_values = strsplit(id_value_str, " ")[[1]]
      rv = c(id_type, id_values)
    }
    return(rv)
  }

  # return true if dblink contains identifier of type <type>
  find_link = function(dblink) {
    rv = FALSE
    parsed_link = parse_link(dblink)
    if(parsed_link[1] == type) {
      rv = TRUE
    }
    return(rv) 
  }

  # reduce dblinks to a (possibly empty) vector of identifiers
  dblinks_to_ids = function(dblinks) {
    ids = NA
    found_link = Find(find_link, dblinks)
    if(length(found_link) == 0) {
      ids = c()
    } else {
      parsed_link = parse_link(found_link)
      ids = parsed_link[2:length(parsed_link)]
    }
    return(ids)
  }

  # keggGet only accepts 10 entries at a time on the server side
  batch_size = 10
  n_gene_ids = length(kegg_gene_ids)
  n_batches = ceiling(n_gene_ids / batch_size)
  for(i in 1:n_batches) {
    start = (i - 1) * batch_size + 1
    end = i * batch_size
    if(end > n_gene_ids) {
      end = n_gene_ids
    }
    batch = kegg_gene_ids[start:end]

    gene_data = keggGet(c(batch))
    gene_ids = mapply(function(gene_datum) {
      return(dblinks_to_ids(gene_datum$DBLINKS))
    }, gene_data)

    link_ids[start:end] = as.array(gene_ids)
  }
  return(link_ids)
}

list_pathways <- function() {
  pathways = keggLink("pathway", "hsa")
  unique_paths = sort(unique(pathways))

  prefix = "path:"
  remove_prefix <- function(str) {
    rv = str
    ind = grep(prefix, str)
    if(!all(ind == 0)) {
      rv = substr(str, ind + nchar(prefix), nchar(str))
    }
    return(rv)
  }

  print(class(unique_paths))
  unique_paths = sapply(unique_paths, remove_prefix)
  return(unique_paths)
}

# TODO upcase
check_id_type = function(id_type) {
  rv = NULL
  if (id_type == 'HGNC') {
    rv = HGNC
  } else if (id_type == 'ENSEMBL') {
    rv = ENSEMBL
  } else if (id_type == 'REFSEQ') {
    rv = REFSEQ
  } else {
    write(paste('--id-type=', id_type, ' is not valid', sep=''), stderr())
    quit(status=3)
  }
  return(rv)
}

main = function() {
  parser = ArgumentParser(description='Download KEGG pathway graphs in the graphml format.')
  parser$add_argument('--id-type', help='Attempt to convert KEGG gene identifiers of form \'hsa:XXXXX\' to this identifier type. Unmapped identifiers will remain as KEGG gene identifiers. Available types: "HGNC", "ENSEMBL", "REFSEQ". Default="ENSEMBL".', default='ENSEMBL')
  parser$add_argument('--outdir', help='Directory to write graphml files to.', required=T)
  args = parser$parse_args()
  id_type = check_id_type(args$id_type)

  any_error = F
  unique_paths = list_pathways()
  for (pathway_id in unique_paths) {
    outfile = file.path(args$outdir, paste(pathway_id, '.graphml', sep=''))
    cat(paste0('Downloading ', pathway_id, '...'), file=stdout())
    if (!file.exists(outfile)) {
      tryCatch(withCallingHandlers({
        graph = get_pathway(pathway_id, args$outdir)
        gene_ids = genes_in_pathway(graph)
        id_map = link_genes_hash(gene_ids, type=id_type)
        graph = rename_nodes(graph, id_map)
        write_graph(graph, outfile, format="graphml")
        write(' done', file=stdout())
      }, error=function(e) {
        write(' failed', file=stdout())
        write(sys.calls(), file=stderr())
        write(message(e), file=stderr())
      }), error=function(e) {
        any_error = T
      })
    } else {
      write(' done', file=stdout())
    }
  }
  if(any_error) {
    quit(status=1)
  } else {
    quit(status=0)
  }
}

main()
