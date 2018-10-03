#!/usr/bin/env Rscript
load_success = TRUE
load_success = load_success && library(KEGGREST, quietly=TRUE, logical.return=TRUE)
load_success = load_success && library(KEGGgraph, quietly=TRUE, logical.return=TRUE)
load_success = load_success && library(graph, quietly=TRUE, logical.return=TRUE)
load_success = load_success && library(igraph, quietly=TRUE, logical.return=TRUE)
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

#' Return a graph object associated with a kegg pathway identifier
#'
#' @param pathway_id e.g. "hsa:05164"
#' @return graph object (see help(graphNEL))
pathway_to_graph <- function(pathway_id) {
  pathway_id_split = strsplit(pathway_id, PATHWAY_DELIM)
  if (pathway_id_split == pathway_id) {
    stop("cannot split %s on %s", pathway_id, PATHWAY_DELIM)
  } 
  pathway_org = pathway_id_split[[1]][1]
  pathway_num = pathway_id_split[[1]][2]
  pathway_file = tempfile()
  kgml = retrieveKGML(pathwayid=pathway_num, organism=pathway_org, destfile=pathway_file)
  kegg.pathway = parseKGML(pathway_file)
  graphnel = KEGGpathway2Graph(kegg.pathway, expandGenes = TRUE)
  graph = graph_from_graphnel(graphnel)
  return(graph)
}

#' Translate KEGG node identifiers using id_map
#'
rename_nodes <- function(G, id_map) {
  nodes = V(G)
  for(node_ind in nodes) {
    node = nodes[node_ind]
    old_name = node$name
    new_name = id_map[[old_name]]
    if(!is.null(new_name)) {
      G = set_vertex_attr(G, 'name', index=old_name, value=new_name)
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

# TODO number of genes off by 1 when compared to pathway_to_graph
genes_in_pathway <- function(pathway_id) {
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

#' Link KEGG gene identifiers e.g. "hsa:5291" to the identifier named in <type>, default ENSEMBL
#'
#' @return ids_map an "environment" used as a hash mapping KEGG gene ids to their ENSEMBL counterparts
#'   or NA if there is no linked identifier
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
#' @return link_ids an array of the same length as kegg_gene_ids with the mapped identifier 
#'   or NA if there is no linked identifier
link_genes_array <- function(kegg_gene_ids, type=ENSEMBL) {
  link_ids = array(dim = length(kegg_gene_ids))

  # return character vector of length 2: <id_type>, <id_value>
  # or an empty vector if link cannot be parsed
  parse_link = function(dblink) {
    rv = c(NULL)
    split_rv = strsplit(dblink, ':', fixed = TRUE)
    if(split_rv == dblink) {
      # then split failed, return NULL
    } else {
      id_type = trimws(split_rv[[1]][1])
      id_value = trimws(split_rv[[1]][2])
      rv = c(id_type, id_value)
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

  # reduce dblinks to a singleton -- an identifier of type <type> or NA
  dblinks_to_id = function(dblinks) {
    id = NA
    found_link = Find(find_link, dblinks)
    if(length(found_link) == 0) {
      id = NA
    } else {
      parsed_link = parse_link(found_link)
      id = parsed_link[2]
    }
    return(id)
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
      return(dblinks_to_id(gene_datum$DBLINKS))
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

main = function() {
  usage = "
Download KEGG pathway graphs in the graphml format.

kegg.R [-h] [-l] [<pathway_id>] [<outfile>]
"
  args <- commandArgs(trailingOnly = TRUE)
  list_flag = FALSE
  help_flag = FALSE

  match_rv = match("-h", args)
  if(!is.na(match_rv)) {
    help_flag = TRUE
  }
  if(help_flag) {
    # then display usage on stdout
    cat(usage, sep="\n")
    quit(status = 0)
  }

  # TODO transform output to e.g. hsa:00010 instead of hsa00010
  match_rv = match("-l", args)
  if(!is.na(match_rv)) {
    list_flag = TRUE
  }
  if(list_flag) {
    # then list the names of all pathways
    unique_paths = list_pathways()
    cat(unique_paths, sep="\n")
    quit(status = 0)
  }

  if(length(args) == 2) {
    # then list all the genes for the given pathway
    pathway_id = args[1] # e.g. "hsa:00010"
    outfile = args[2]
    gene_ids = genes_in_pathway(pathway_id)
    id_map = link_genes_hash(gene_ids, type=HGNC)
    graph = pathway_to_graph(pathway_id)
    graph = rename_nodes(graph, id_map)
    write_graph(graph, outfile, format="graphml")
    #write_abc(graph, outfile, id_map = id_map)
  } else {
    write(usage, stderr()) 
    quit(status = 1)
  }
}

main()
