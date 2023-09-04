#include "SelectionBase.h"

#include <iostream>
#include "Functions.h"

SelectionBase::SelectionBase(std::string fSelectionName_) {
  fSelectionName = fSelectionName_;
  nPassedEvents = 0;
}

void SelectionBase::Setup(TTree* Tree_, bool Create_) {
  SetupTree(Tree_, Create_);
  DefineConstants();
}

void SelectionBase::ApplySelection(AnalysisEvent* Event) {
  Reset();

  Selected = Selection(Event);
  MC_Signal = DefineSignal(Event);
  EventCategory = categorize_event(Event, TrueFV);
  
  ComputeObservables(Event);
  
  if (Selected) {
    nPassedEvents++;
  }
}

void SelectionBase::Print() {
  std::cout << fSelectionName << " has " << nPassedEvents << " events which passed" << std::endl;
}

void SelectionBase::SetupTree(TTree* Tree_, bool Create_) {
  Tree = Tree_;
  Create = Create_;

  std::string BranchName;

  BranchName = "Selected";
  SetBranch(&Selected,BranchName,kBool);

  BranchName = "MC_Signal";
  SetBranch(&MC_Signal,BranchName,kBool);

  BranchName = fSelectionName+"_Category";
  SetBranch(&EventCategory,"EventCategory",kInteger);
  
  DefineBranches();
}

void SelectionBase::SetBranch(void* Variable, std::string VariableName, VarType VariableType) {
  SaveVariablePointer(Variable,VariableType);

  VariableName = fSelectionName+"_"+VariableName;
  std::string Leaflist = VariableName;

  switch (VariableType) {
  case kBool:
    Leaflist += "/O";
    break;
  case kDouble:
    Leaflist += "/D";
    break;
  case kFloat:
    Leaflist += "/F";
    break;
  case kInteger:
    Leaflist += "/I";
    break;
  default:
    std::cerr << "Unexpected variable type:" << VariableType << std::endl;
    throw;
  }

  set_output_branch_address(*Tree,VariableName,Variable,Create,Leaflist);
}

void SelectionBase::SaveVariablePointer(void* Variable, VarType VariableType) {
  switch (VariableType) {
  case kBool:
    Pointer_Bool.push_back((bool*)Variable);
    break;
  case kDouble:
    Pointer_Double.push_back((double*)Variable);
    break;
  case kFloat:
    Pointer_Float.push_back((float*)Variable);
    break;
  case kInteger:
    Pointer_Integer.push_back((int*)Variable);
    break;
  default:
    std::cerr << "Unexpected variable type:" << VariableType << std::endl;
    throw;
  }
}

void SelectionBase::Reset() {
  for (size_t i=0;i<Pointer_Bool.size();i++) {
    *(Pointer_Bool[i]) = false;
  }
  for (size_t i=0;i<Pointer_Double.size();i++) {
    *(Pointer_Double[i]) = BOGUS;
  }
  for (size_t i=0;i<Pointer_Float.size();i++) {
    *(Pointer_Float[i]) = BOGUS;
  }
  for (size_t i=0;i<Pointer_Integer.size();i++) {
    *(Pointer_Integer[i]) = BOGUS_INDEX;
  }
}
