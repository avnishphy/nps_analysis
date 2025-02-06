#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <iostream>
#include <fstream>

void ListKeysBranchesLeaves(const char* filename, const char* outputfile) {
    // Open the ROOT file
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    // Open output file
    std::ofstream outfile(outputfile);
    if (!outfile.is_open()) {
        std::cerr << "Error: Cannot open output file " << outputfile << std::endl;
        return;
    }

    outfile << "ROOT File: " << filename << std::endl << std::endl;

    // List all keys in the ROOT file
    outfile << "Keys in the file:" << std::endl;
    TIter nextkey(file->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)nextkey())) {
        outfile << "Key Name: " << key->GetName() << ", Class: " << key->GetClassName() << std::endl;

        // Check if the key is a tree
        if (strcmp(key->GetClassName(), "TTree") == 0) {
            TTree* tree = (TTree*)file->Get(key->GetName());
            if (tree) {
                outfile << "  Branches in tree " << key->GetName() << ":" << std::endl;

                // Iterate over branches
                TIter nextbranch(tree->GetListOfBranches());
                TBranch* branch;
                while ((branch = (TBranch*)nextbranch())) {
                    outfile << "    Branch Name: " << branch->GetName() << std::endl;

                    // Iterate over leaves in the branch
                    TIter nextleaf(branch->GetListOfLeaves());
                    TLeaf* leaf;
                    while ((leaf = (TLeaf*)nextleaf())) {
                        outfile << "      Leaf Name: " << leaf->GetName() << ", Type: " << leaf->GetTypeName() << std::endl;
                    }
                }
            }
        }
    }

    // Close files
    outfile.close();
    file->Close();

    std::cout << "Keys, branches, and leaves saved to: " << outputfile << std::endl;
}
