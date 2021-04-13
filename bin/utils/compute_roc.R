

compute_roc <- function(scores, labels, outdir){

    # create data
    dat = mmdata(scores, labels)

    # Calculate ROC and Precision-Recall curves for multiple models
    mscurves <- evalmod(dat)

    # Plot and save ROC and Precision-Recall curves 
    jpeg(file=file.path(outdir, "curve.jpg"))
    ggplot2::autoplot(mscurves)
    dev.off()

    # Retrieve a dataframe of AUCs
    aucs.df <- auc(mscurves)
    roc = aucs.df[which(aucs.df$curvetypes=="ROC"),"aucs"]
    pr = aucs.df[which(aucs.df$curvetypes=="PRC"),"aucs"]
    writeLines(as.character(roc), file.path(outdir, "roc.txt"))
    writeLines(as.character(pr), file.path(outdir, "pr.txt"))
}

