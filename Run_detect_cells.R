#/usr/bin/env Rscript

require(DropletUtils)
require(Seurat)
require(optparse)


option_list <- list(
  make_option(c("-p", "--path_file"), action="store",
              help="path file for 10X h5 files", type="character"),
  make_option(c("-f", "--fdr"), action="store", default=0.01,
              help="False Discovery Rate Cutoff", type="double"),
  make_option(c("-i", "--iters"), action="store", default=10000,
              help="Number of iterations for DetectCells()",type="integer"),
  make_option(c("-o", "--out_path"), action="store", default=".",
              help="Directory to save the reports and matrices",type="character"),
  make_option(c("-s", "--save_matrix"), action="store", default = TRUE, type = "logical",
              help="Save the FDR filtered matrix (in 10X format)"),
  make_option(c("-r", "--write_report"), action="store", default = TRUE, type="logical",
              help="Write an HTML report")

)

args = parse_args(OptionParser(option_list=option_list))




paths = readLines(args$path_file)


RunDetectCells = function(path, out_path=args$out_path, iters=args$iters, fdr=args$fdr, save_matrix=args$save_matrix, write_report=args$write_report){
  
  try_names = strsplit(path, "/")[[1]]
  if(tail(try_names, n=1) == "outs"){
    sample = tail(try_names, n=2)[1]
  } else{
    sample = tail(try_names, n=1)
  }
   
  
  print(paste0("Reading h5 for sample: ", sample))
  h5 = Read10X_h5(paste0(path,"raw_feature_bc_matrix.h5"))
  
  print(paste0("Running emtyDrops for sample: ", sample))
  set.seed(100)
  e.out <- emptyDrops(h5, niters=iters)
  is.cell <- e.out$FDR <= fdr 
  is.cell[is.na(is.cell)] = FALSE
  fil_matrix = h5[,is.cell] 
  
  if(save_matrix){
    print(paste0("writing filtered matrix file for: ",sample))
    write10xCounts(x = fil_matrix, path = paste0(out_path,"/",sample,"/"),type = "sparse", overwrite = TRUE)
  }else{}
  
  if(write_report){
    print(paste0("writing markdown report for: ",sample))
    if(!dir.exists(paste0(out_path,"/reports"))){
      dir.create(paste0(out_path,"/reports"))
    }
    rmarkdown::render('/volumes/lab/users/aschalck/TENX/analysis/BCR/BCR_4/detect_cells_report.Rmd',  # file 2
                      output_file =  paste0(sample, '_detect_cells_report.html'), 
                      output_dir = paste0(out_path,"/reports")) 
  }else{}
  
}

lapply(paths, function(x) RunDetectCells(path=x))


 


