

library(tidyverse)
library(fuzzyjoin)
library(pmultinom)
library(parallel)


#clonotype = rbind(clonotypes_list$BCR20_1,clonotypes_list$BCR20_2)



#merge clonotypes function
merge_clonotypes = function(clonotype){
  
  #all clonotypes with alpha and beta
  clonotype %>%
    filter(str_detect(cdr3s_nt, pattern = "TRA:") & str_detect(cdr3s_nt, pattern = "TRB:")) ->
    complete
  
  #just alpha clonotypes
  clonotype %>%
    filter(str_detect(cdr3s_nt, pattern = "TRA:") & !str_detect(cdr3s_nt, pattern = "TRB:")) %>% 
    mutate(cdr3s_nt = as.character(cdr3s_nt))->
    alpha_only
  
  #just beta clonotypes
  clonotype %>%
    filter(!str_detect(cdr3s_nt, pattern = "TRA:") & str_detect(cdr3s_nt, pattern = "TRB:")) %>% 
    mutate(cdr3s_nt = as.character(cdr3s_nt)) ->
    beta_only
  
  #annotate just beta clonotpes to matching complete
  complete.beta=fuzzy_join(complete,beta_only, by = "cdr3s_nt", match_fun = str_detect, mode="left")
  names(complete.beta) = c(names(complete), paste0(names(beta_only), ".beta"))
  
  #annotate just alpha clonotpes to matching complete
  complete.alpha=fuzzy_join(complete,alpha_only, by = "cdr3s_nt", match_fun = str_detect, mode="left")
  names(complete.alpha) = c(names(complete), paste0(names(alpha_only), ".alpha"))
  
  #merge the beta and alpha merges into one with both beta and alpha singlets per complete
  #merged = full_join(complete.beta, complete.alpha)
  
  #merge the incomplete betas with just the top complete clonotype
  complete.beta %>% 
    filter(!is.na(clonotype_id.beta)) %>% 
    group_by(clonotype_id.beta) %>% 
    top_n(n=1, wt=frequency) %>%
    top_n(n=-1, wt=clonotype_id) %>% 
    ungroup() ->
    complete.beta.only
  
  #merge the incomplete alphas with just the top complete clonotype
  complete.alpha %>% 
    filter(!is.na(clonotype_id.alpha)) %>% 
    group_by(clonotype_id.alpha) %>% 
    top_n(n=1, wt=frequency) %>%
    top_n(n=-1, wt=clonotype_id) %>% 
    ungroup() ->
    complete.alpha.only
  
  #add the merged columns back to the complete table
  #label the rows as complete clonotypes
  #label the rows as merged or not
  complete %>%
    left_join(complete.beta.only) %>% 
    left_join(complete.alpha.only) %>%
    mutate(merge="unmerged", type="complete") %>% 
    mutate(merge=ifelse(!is.na(clonotype_id.alpha),"alpha",merge)) %>% 
    mutate(merge=ifelse(!is.na(clonotype_id.beta),"beta",merge)) %>% 
    mutate(merge=ifelse((!is.na(clonotype_id.alpha) & !is.na(clonotype_id.beta)),"both",merge)) ->
    complete.labeled
  
  #find all the clonotypes associated with a compelte clonotype
  complete.clonotypes = c(as.character(complete.labeled$clonotype_id), as.character(complete.labeled$clonotype_id.beta), as.character(complete.labeled$clonotype_id.alpha))
  complete.clonotypes = complete.clonotypes[!is.na(complete.clonotypes)]
  
  #find alpha onlys that don't get merged
  alpha_only %>% 
    filter(!clonotype_id %in% complete.clonotypes) %>% 
    mutate(merge="unmerged", type="alpha_only") ->
    single.alpha
  
  #find beta onlys that don't get merged
  beta_only %>% 
    filter(!clonotype_id %in% complete.clonotypes) %>% 
    mutate(merge="unmerged", type="beta_only") ->
    single.beta
  
  #add the unmerged betas and alphas back
  complete.labeled %>% 
    full_join(single.beta) %>% 
    full_join(single.alpha) %>%
    mutate(frequency.beta = if_else(is.na(frequency.beta),0,as.numeric(frequency.beta))) %>% 
    mutate(frequency.alpha = if_else(is.na(frequency.alpha),0,as.numeric(frequency.alpha))) %>% 
    mutate(proportion.beta = if_else(is.na(proportion.beta),0,as.numeric(proportion.beta))) %>% 
    mutate(proportion.alpha = if_else(is.na(proportion.alpha),0,as.numeric(proportion.alpha))) -> 
    clonotype_merged_unsummed
    
  #sum the new frequencies and proportions beta
  clonotype_merged_unsummed %>% 
    filter(!is.na(clonotype_id.beta)) %>% 
    filter(!duplicated(clonotype_id.beta)) %>% 
    group_by(clonotype_id) %>%
    summarise(frequency.beta.new = sum(frequency.beta), proportion.beta.new = sum(proportion.beta))->
    clonotype_merged_summed.beta
  
  #sum the new frequencies and proportions alpha
  clonotype_merged_unsummed %>% 
    filter(!is.na(clonotype_id.alpha)) %>% 
    filter(!duplicated(clonotype_id.alpha)) %>% 
    group_by(clonotype_id) %>%
    summarise(frequency.alpha.new = sum(frequency.alpha), proportion.alpha.new = sum(proportion.alpha)) ->
    clonotype_merged_summed.alpha
  
  #join back the sums and final sum
  clonotype_merged_unsummed %>% 
    left_join(clonotype_merged_summed.beta) %>% 
    left_join(clonotype_merged_summed.alpha) %>% 
    mutate(frequency.beta.new = if_else(is.na(frequency.beta.new),0,as.numeric(frequency.beta.new))) %>% 
    mutate(frequency.alpha.new = if_else(is.na(frequency.alpha.new),0,as.numeric(frequency.alpha.new))) %>% 
    mutate(proportion.beta.new = if_else(is.na(proportion.beta.new),0,as.numeric(proportion.beta.new))) %>% 
    mutate(proportion.alpha.new = if_else(is.na(proportion.alpha.new),0,as.numeric(proportion.alpha.new))) %>% 
    mutate(frequency.new = frequency + frequency.beta.new + frequency.alpha.new) %>% 
    mutate(proportion.new = proportion + proportion.beta.new + proportion.alpha.new) ->
    clonotype_merged
  
  return(clonotype_merged)
}

#annontate expanded clonotypes function
annotate_expanded = function(clonotype_merged, fwer){
  
  #number of clones and cells
  nclones = length(unique(clonotype_merged$clonotype_id))
  clonotype_merged %>% 
    filter(!duplicated(clonotype_id)) %>%
    pull(frequency.new) %>% 
    sum() -> ncells

  #null distribution where they are all the same
  theoretical.dist <- rep(1/nclones,nclones)
  
  #pvalue cutoff function
  cutoff <- Vectorize(function(prob, target.pval, size)
  {
    pvals <- 1 - pbinom(0:size, size, prob)
    min((0:size)[pvals <= target.pval])
  })
  
  # p values we want to try (just 100 values between 0 and .05)
  xvals <- seq.int(0, fwer, length.out=100)
  
  # FWER for each of these values
  yvals <- lapply(xvals, function(p) {1 - pmultinom(upper=cutoff(theoretical.dist, p, ncells), size=ncells, probs=theoretical.dist, method="exact"); print(p)})
  
  #plot(xvals, yvals, type='l',
  #     main="Controlling family-wise error rate", xlab="binomial p-value", ylab="FWER")
  
  #choose pvalue less than fwer
  maxp = max(xvals[yvals <= fwer])
  
  #neg biom pvalues per cell
  pexpand = lapply(clonotype_merged$frequency.new, function(x) 1-pbinom(q = x,size = ncells, prob = 1/nclones))
  
  #true false annotation
  pexpand_TF = if_else(pexpand <= maxp, "expanded", "unexpanded")
  
  #annotate the original df
  clonotype_merged$p_expand = pexpand
  clonotype_merged$expand = pexpand_TF
  
  
  #note the top expanded clonotype(s)
  clonotype_merged %>% 
    mutate(top_clones = if_else(frequency.new==max(frequency.new),"yes","no")) ->
    clonotype_merged
  
  return(clonotype_merged)
}


#aca = rbind(contig_anno_list$BCR20_1, contig_anno_list$BCR20_2)
#cme = merge_clonotypes(clonotype)


#annontate the cells function
annotate_cells = function(aca,cme){
  aca %>% 
    select(barcode,is_cell,raw_clonotype_id) %>% 
    filter(raw_clonotype_id!="None") %>% 
    unique() -> 
    aca_fil
  
  cme %>% 
    select(clonotype_id,clonotype_id.beta) %>% 
    unique() ->
    cme.beta
  
  cme %>% 
    select(clonotype_id,clonotype_id.alpha) %>% 
    unique() ->
    cme.alpha 
  
  left_join(aca_fil, cme.beta, by=c("raw_clonotype_id"="clonotype_id.beta")) %>%
    left_join(cme.alpha, by=c("raw_clonotype_id"="clonotype_id.alpha"), suffix=c(".p.beta",".p.alpha")) %>%
    mutate(new.clonotype_id = raw_clonotype_id) %>% 
    mutate(new.clonotype_id = if_else(!is.na(clonotype_id.p.beta),as.character(clonotype_id.p.beta),as.character(new.clonotype_id))) %>% 
    mutate(new.clonotype_id = if_else(!is.na(clonotype_id.p.alpha),as.character(clonotype_id.p.alpha),as.character(new.clonotype_id))) %>% 
    mutate(cell_merge = "none") %>% 
    mutate(cell_merge = if_else(!is.na(clonotype_id.p.beta),"beta",cell_merge)) %>% 
    mutate(cell_merge = if_else(!is.na(clonotype_id.p.alpha),"alpha",cell_merge))  %>% 
    select(-clonotype_id.p.beta,-clonotype_id.p.alpha) ->
    aca.beta.alpha
  
  
  left_join(aca.beta.alpha, cme, by=c("new.clonotype_id"="clonotype_id")) %>% 
    select(-contains("alpha"),-contains("beta")) %>% 
    unique() ->
    #mutate(barcode.new = paste0(name,"___",str_sub(barcode,start = 1,end = 16))) %>% 
    #mutate(raw_clonotype_id=paste0(raw_clonotype_id,"-",name)) %>% 
    #mutate(new.clonotype_id=paste0(new.clonotype_id,"-",name)) 
    cell.clone.merge
  
  return(cell.clone.merge)
  
}


#test = annotate_cells(aca, cme)







































