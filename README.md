# 2pCorrelation
To run this:
# Step1: 
Upload this to your STAR account
# Step2: 
Run the following commands

  root -l  
  .L AnaCorrFourier.C+  
  AnaCorrFourier* agp = new AnaCorrFourier()  
  agp->from = 0  
  agp->to = 2535                      // Can be any number of events you want to run it on. 
  agp->filelist = "input2.txt"        // txt file containing your data file location and name. 
  agp->outfile  = "Analysis.root"  
  agp->SetTotal(2535)  
  agp->SetZbinSize(1)  
  agp->run()  
  
# if the previous snippet runs write a bash code for job submission
