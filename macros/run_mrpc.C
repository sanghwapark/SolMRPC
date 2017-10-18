void run_mrpc(string fin, string fout, int ngen)
{

  gSystem->Load("libsolmrpc");

  SolMRPCdigit* sd = new SolMRPCdigit("SolMRPC",fin);
  sd->SetOutFileName(fout);
  sd->SetSourceMode(0);
  sd->Init();
  //  sd->DoSingle(0);                                                                                                                
  //  sd->SetQthreshold(40.e-3);                                                                                                      
  sd->SetQthreshold(15.e-3); //15fC                                                                                                   
  sd->process(0, ngen); // 1: # of events you want to process, 2: # of total generated events                                         
  sd->End();

  delete sd;
}
