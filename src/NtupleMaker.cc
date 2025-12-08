// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      NtupleMaker
// 
// Description: Dump mumugamma decays
//
// Author:  Zhen Hu
//
/*
    Version: Monte-Carlo 2018 Ntuple(Bdecay)
    Modified by: Shiyang CHEN
*/

// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/Common/interface/View.h>
//For kinematic fit:
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include "TH2F.h"

#define VERBOSE// switch between verbose output and concise output
#ifdef VERBOSE
    #define vout cout
#else
    #define vout {}/##/
#endif
#define SQUARE(x) ((x) * (x))

using namespace std;
using namespace edm;
using namespace reco;
using namespace muon;
using namespace trigger;
// mass constant
const float muonMass = 0.1056583;
float muonSigma = 0.000001;// KinematicParticleFactory#@$%.particle(,,,,float&) left value ref cannot bind constants
const float JpsiMass = 3.096900;
const float psi2SMass = 3.686097;
const float JpsiSigma = 0.00004;
const float psi2SSigma = 0.00004;
//
// class declaration
//
class NtupleMaker: public edm::EDAnalyzer {
    public:
        explicit NtupleMaker(const edm::ParameterSet &);
        ~NtupleMaker();
        static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

    private:
        UInt_t getTriggerBits(const edm::Event &);
        bool triggerMatch(TLorentzVector mu1p4, TLorentzVector mu2p4, TLorentzVector mu3p4, TLorentzVector mu4p4);
        bool findTrigger(edm::Handle<edm::TriggerResults>& hltR, vector < string >& triggersToCheck, vector < string >& triggerNameFound);
        void analyzeTrigger(edm::Handle<edm::TriggerResults>& hltR, edm::Handle<trigger::TriggerEvent>& hltE, const string& triggerName);
        bool triggerDecision(edm::Handle<edm::TriggerResults>& hltR, int iTrigger);

        void setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles, const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag>& pv);
        void fourMuonFit(edm::Handle< edm::View<pat::Muon> > muons, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV, edm::Handle<reco::VertexCollection>);
        void calCtau(RefCountedKinematicVertex& decayVrtx, RefCountedKinematicParticle& kinePart, Vertex& bs, Double_t& ctau, Double_t& ctauErr, Double_t& LxyPV, Double_t& sigLxy);// dedicated to RECO level
        double calCtau(const reco::Candidate& p, const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag>& pV);// dedicated to GEN level
        Double_t calD(RefCountedKinematicVertex decayVrtx1, RefCountedKinematicVertex decayVrtx2);

        virtual void beginJob();
        virtual void analyze(const edm::Event &, const edm::EventSetup &);
        virtual void endJob();

        virtual void beginRun(edm::Run &, edm::EventSetup const&);
        virtual void endRun(edm::Run const &, edm::EventSetup const &);
        virtual void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);
        virtual void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);
		// ----------member data ---------------------------
        edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_Label;
		// edm::EDGetTokenT<pat::CompositeCandidateCollection> conversion_Label;
        edm::EDGetTokenT<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag>> primaryVtxToken_;
		edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
		edm::EDGetTokenT<reco::BeamSpot> bs_Label;
		edm::EDGetTokenT<edm::View<pat::Muon>> muon_Label;
		edm::EDGetTokenT<edm::TriggerResults> triggerResultsTok_;
        edm::EDGetTokenT<trigger::TriggerEvent> triggerEventTok_;
        //edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

        // // tree data arranged in struct
        // // -------------------------
        // // GEN level information
        // // -------------------------
        // vector<struct GENmu> GenMuonList;
        // struct GENdimu GenJpsi, Genpsi2S;
        // struct GENevt GenEvent;
        // // -------------------------
        // // RECO level information
        // // -------------------------
        // vector<struct RECOmu> RecoMuonList;
        // vector<struct RECOdimu> RecoJpsiList, Recopsi2SList;
        // vector<struct RECOevt> RecoEventList;

        // GEN level muon
        vector<Double_t> GEmu_pt, GEmu_eta, GEmu_phi, GEmu_mass;
        vector<int> GEmu_recoId;
        // GEN level J/psi
        Double_t GEJpsi_pt, GEJpsi_eta, GEJpsi_y, GEJpsi_phi, GEJpsi_mass, GEJpsi_ctau;
        int GEJpsi_recoId, GEJpsi_muId1, GEJpsi_muId2;
        bool GEJpsi_matchGEN, GEJpsi_passID, GEJpsi_passVtx;
        // GEN level psi(2S)
        Double_t GEpsi2S_pt, GEpsi2S_eta, GEpsi2S_y, GEpsi2S_phi, GEpsi2S_mass, GEpsi2S_ctau;
        int GEpsi2S_recoId, GEpsi2S_muId1, GEpsi2S_muId2;
        bool GEpsi2S_matchGEN, GEpsi2S_passID, GEpsi2S_passVtx;
        // GEN level event
        vector<double> GEevt_vtx;
        Double_t GEevt_fourMuMass;
        bool GEevt_valid, GEevt_passAcc, GEevt_passHLT, GEevt_passAll;
        // RECO level muon
        vector<Double_t> REmu_pt, REmu_eta, REmu_phi, REmu_mass;
        vector<bool> REmu_isSoft, REmu_passCut;
        // RECO level J/psi
        vector<Double_t> REJpsi_pt, REJpsi_eta, REJpsi_y, REJpsi_phi, REJpsi_mass, REJpsi_massErr, REJpsi_vtxProb, REJpsi_cstrVtxProb,
        REJpsi_ctau, REJpsi_ctauErr, REJpsi_LxyPV, REJpsi_sigLxy;
        vector<int> REJpsi_muId1, REJpsi_muId2;
        vector<bool> REJpsi_passCut;
        // RECO level psi(2S)
        vector<Double_t> REpsi2S_pt, REpsi2S_eta, REpsi2S_y, REpsi2S_phi, REpsi2S_mass, REpsi2S_massErr, REpsi2S_vtxProb, REpsi2S_cstrVtxProb,
        REpsi2S_ctau, REpsi2S_ctauErr, REpsi2S_LxyPV, REpsi2S_sigLxy;
        vector<int> REpsi2S_muId1, REpsi2S_muId2;
        vector<bool> REpsi2S_passCut;
        // RECO level event
        vector<double> REevt_x, REevt_y, REevt_z;
        vector<Double_t> REevt_fourMuMass, REevt_massChisq, REevt_vtxProb, REevt_d;
        vector<int> REevt_JpsiId, REevt_psi2SId;
        vector<bool> REevt_passHLT, REevt_matchTrg, REevt_fourMuFit, REevt_twoDimuFit;

        int Total_events_analyzed;
        int Total_events_analyzed_triggered;
        int Total_events_dimuon;
        int Total_events_dimuon_trg_matched;
        int Total_events_Fourmuon;
        // some counters for cut flow
        int nevent_soft4mu;
        int nevent_4mu_pt;
        int nevent_4mu_eta;
        int nevent_dimu_vtx;
        int nevent_fourmu; 
        int nevent_jpsi_pT;
		ULong64_t run;
		ULong64_t lumi;
        ULong64_t event;
		UInt_t trigger;
 
        string hltName_;
        string   triggerName_;
        edm::TriggerNames triggerNames;
        vector<string> triggerList;
        int runNumber;
  		// vector<reco::MuonCollection::const_iterator> allL1TrigMuons;
  		// vector<reco::MuonCollection::const_iterator> allL2TrigMuons;
  		// vector<reco::MuonCollection::const_iterator> allL3TrigMuons;
  		// vector<GlobalVector> allMuL1TriggerVectors;
  		// vector<GlobalVector> allMuL2TriggerVectors;
  		// vector<GlobalVector> allMuL3TriggerVectors_lowEff;
  		// vector<GlobalVector> allMuL3TriggerVectors_highEff;
  		vector<GlobalVector> muonTriggeredVecs;
        vector<GlobalVector> JpsiTriggeredVecs;
        HLTConfigProvider hltConfig_;
		TTree *onia_tree;

        edm::EDGetTokenT<reco::GenParticleCollection> prunedgenParticlesSrc_;
        //edm::EDGetTokenT<pat::PackedGenParticle> packedgenParticlesSrc_;
        // Counters
};
//
// constructors and destructor
//
NtupleMaker::NtupleMaker(const edm::ParameterSet & iConfig):
    dimuon_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("dimuons"))),
    // conversion_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("conversions"))),
    primaryVtxToken_(mayConsume<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag>>(iConfig.getParameter<edm::InputTag>("primaryVertex"))),
    primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
    bs_Label(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("offlineBeamSpot"))),
    muon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter< edm::InputTag>("muons"))),
    triggerResultsTok_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
    triggerEventTok_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("TriggerSummaryAOD"))),
    triggerList(iConfig.getUntrackedParameter<vector<string>>("triggerList")) {
	
    edm::Service < TFileService > fs;
	onia_tree = fs->make < TTree > ("oniaTree", "Tree of MuMuGamma");
    Total_events_analyzed = 0;
    nevent_soft4mu = 0;
    nevent_4mu_pt=0;
    nevent_4mu_eta=0;
    Total_events_analyzed_triggered = 0 ;
    Total_events_dimuon= 0;
    Total_events_dimuon_trg_matched = 0;
    Total_events_Fourmuon = 0;

    runNumber = -99;
    hltName_ = "HLT";
    triggerName_ = "@"; // "@" means: analyze all triggers in config

    onia_tree->Branch("run", &run, "run/l");
    onia_tree->Branch("lumi", &lumi, "lumi/l");
    onia_tree->Branch("event", &event, "event/l");
    onia_tree->Branch("trigger", &trigger, "trigger/I");
    // GEN level
    onia_tree->Branch("GEmu_pt", &GEmu_pt);
    onia_tree->Branch("GEmu_eta", &GEmu_eta);
    onia_tree->Branch("GEmu_phi", &GEmu_phi);
    onia_tree->Branch("GEmu_mass", &GEmu_mass);
    onia_tree->Branch("GEmu_recoId", &GEmu_recoId);
    // GEN level J/psi
    onia_tree->Branch("GEJpsi_pt", &GEJpsi_pt, "GEJpsi_pt/D");
    onia_tree->Branch("GEJpsi_eta", &GEJpsi_eta, "GEJpsi_eta/D");
    onia_tree->Branch("GEJpsi_y", &GEJpsi_y, "GEJpsi_y/D");
    onia_tree->Branch("GEJpsi_phi", &GEJpsi_phi, "GEJpsi_phi/D");
    onia_tree->Branch("GEJpsi_mass", &GEJpsi_mass, "GEJpsi_mass/D");
    onia_tree->Branch("GEJpsi_ctau", &GEJpsi_ctau, "GEJpsi_ctau/D");
    onia_tree->Branch("GEJpsi_muId1", &GEJpsi_muId1, "GEJpsi_muId1/I");
    onia_tree->Branch("GEJpsi_muId2", &GEJpsi_muId2, "GEJpsi_muId2/I");
    onia_tree->Branch("GEJpsi_recoId", &GEJpsi_recoId, "GEJpsi_recoId/I");
    onia_tree->Branch("GEJpsi_matchGEN", &GEJpsi_matchGEN, "GEJpsi_matchGEN/b");
    onia_tree->Branch("GEJpsi_passID", &GEJpsi_passID, "GEJpsi_passID/b");
    onia_tree->Branch("GEJpsi_passVtx", &GEJpsi_passVtx, "GEJpsi_passVtx/b");
    // GEN level psi(2S)
    onia_tree->Branch("GEpsi2S_pt", &GEpsi2S_pt, "GEpsi2S_pt/D");
    onia_tree->Branch("GEpsi2S_eta", &GEpsi2S_eta, "GEpsi2S_eta/D");
    onia_tree->Branch("GEpsi2S_y", &GEpsi2S_y, "GEpsi2S_y/D");
    onia_tree->Branch("GEpsi2S_phi", &GEpsi2S_phi, "GEpsi2S_phi/D");
    onia_tree->Branch("GEpsi2S_mass", &GEpsi2S_mass, "GEpsi2S_mass/D");
    onia_tree->Branch("GEpsi2S_ctau", &GEpsi2S_ctau, "GEpsi2S_ctau/D");
    onia_tree->Branch("GEpsi2S_muId1", &GEpsi2S_muId1, "GEpsi2S_muId1/I");
    onia_tree->Branch("GEpsi2S_muId2", &GEpsi2S_muId2, "GEpsi2S_muId2/I");
    onia_tree->Branch("GEpsi2S_recoId", &GEpsi2S_recoId, "GEpsi2S_recoId/I");
    onia_tree->Branch("GEpsi2S_matchGEN", &GEpsi2S_matchGEN, "GEpsi2S_matchGEN/b");
    onia_tree->Branch("GEpsi2S_passID", &GEpsi2S_passID, "GEpsi2S_passID/b");
    onia_tree->Branch("GEpsi2S_passVtx", &GEpsi2S_passVtx, "GEpsi2S_passVtx/b");
    // GEN level event
    onia_tree->Branch("GEevt_vtx", &GEevt_vtx);
    onia_tree->Branch("GEevt_fourMuMass", &GEevt_fourMuMass, "GEevt_fourMuMass/D");
    onia_tree->Branch("GEevt_valid", &GEevt_valid, "GEevt_valid/b");
    onia_tree->Branch("GEevt_passAcc", &GEevt_passAcc, "GEevt_passAcc/b");
    onia_tree->Branch("GEevt_passHLT", &GEevt_passHLT, "GEevt_passHLT/b");
    onia_tree->Branch("GEevt_passAll", &GEevt_passAll, "GEevt_passAll/b");
    // RECO level muon
    onia_tree->Branch("REmu_pt", &REmu_pt);
    onia_tree->Branch("REmu_eta", &REmu_eta);
    onia_tree->Branch("REmu_phi", &REmu_phi);
    onia_tree->Branch("REmu_mass", &REmu_mass);
    onia_tree->Branch("REmu_isSoft", &REmu_isSoft);
    onia_tree->Branch("REmu_passCut", &REmu_passCut);
    // RECO level J/psi
    onia_tree->Branch("REJpsi_pt", &REJpsi_pt);
    onia_tree->Branch("REJpsi_eta", &REJpsi_eta);
    onia_tree->Branch("REJpsi_y", &REJpsi_y);
    onia_tree->Branch("REJpsi_phi", &REJpsi_phi);
    onia_tree->Branch("REJpsi_mass", &REJpsi_mass);
    onia_tree->Branch("REJpsi_massErr", &REJpsi_massErr);
    onia_tree->Branch("REJpsi_vtxProb", &REJpsi_vtxProb);
    onia_tree->Branch("REJpsi_cstrVtxProb", &REJpsi_cstrVtxProb);
    onia_tree->Branch("REJpsi_ctau", &REJpsi_ctau);
    onia_tree->Branch("REJpsi_ctauErr", &REJpsi_ctauErr);
    onia_tree->Branch("REJpsi_LxyPV", &REJpsi_LxyPV);
    onia_tree->Branch("REJpsi_sigLxy", &REJpsi_sigLxy);
    onia_tree->Branch("REJpsi_muId1", &REJpsi_muId1);
    onia_tree->Branch("REJpsi_muId2", &REJpsi_muId2);
    onia_tree->Branch("REJpsi_passCut", &REJpsi_passCut);
    // RECO level psi(2S)
    onia_tree->Branch("REpsi2S_pt", &REpsi2S_pt);
    onia_tree->Branch("REpsi2S_eta", &REpsi2S_eta);
    onia_tree->Branch("REpsi2S_y", &REpsi2S_y);
    onia_tree->Branch("REpsi2S_phi", &REpsi2S_phi);
    onia_tree->Branch("REpsi2S_mass", &REpsi2S_mass);
    onia_tree->Branch("REpsi2S_massErr", &REpsi2S_massErr);
    onia_tree->Branch("REpsi2S_vtxProb", &REpsi2S_vtxProb);
    onia_tree->Branch("REpsi2S_cstrVtxProb", &REpsi2S_cstrVtxProb);
    onia_tree->Branch("REpsi2S_ctau", &REpsi2S_ctau);
    onia_tree->Branch("REpsi2S_ctauErr", &REpsi2S_ctauErr);
    onia_tree->Branch("REpsi2S_LxyPV", &REpsi2S_LxyPV);
    onia_tree->Branch("REpsi2S_sigLxy", &REpsi2S_sigLxy);
    onia_tree->Branch("REpsi2S_muId1", &REpsi2S_muId1);
    onia_tree->Branch("REpsi2S_muId2", &REpsi2S_muId2);
    onia_tree->Branch("REpsi2S_passCut", &REpsi2S_passCut);
    // RECO level event
    onia_tree->Branch("REevt_x", &REevt_x);
    onia_tree->Branch("REevt_y", &REevt_y);
    onia_tree->Branch("REevt_z", &REevt_z);
    onia_tree->Branch("REevt_fourMuMass", &REevt_fourMuMass);
    onia_tree->Branch("REevt_massChisq", &REevt_massChisq);
    onia_tree->Branch("REevt_vtxProb", &REevt_vtxProb);
    onia_tree->Branch("REevt_d", &REevt_d);
    onia_tree->Branch("REevt_JpsiId", &REevt_JpsiId);
    onia_tree->Branch("REevt_psi2SId", &REevt_psi2SId);
    onia_tree->Branch("REevt_passHLT", &REevt_passHLT);
    onia_tree->Branch("REevt_matchTrg", &REevt_matchTrg);
    onia_tree->Branch("REevt_fourMuFit", &REevt_fourMuFit);
    onia_tree->Branch("REevt_twoDimuFit", &REevt_twoDimuFit);

	prunedgenParticlesSrc_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
	//packedgenParticlesSrc_ = consumes<pat::PackedGenParticle>((edm::InputTag)"packedGenParticles");
    //prunedgenParticlesSrc_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("prunedgenParticlesSrc")));
    //packedgenParticlesSrc_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("packedgenParticlesSrc")));
}

NtupleMaker::~NtupleMaker() {}
//
// member functions
//
/*
    Grab Trigger information. Save it in variable trigger, trigger is an int between 0 and 127, in binary it is:
    (pass 2)(pass 1)(pass 0)
    ex. 7 = pass 0, 1 and 2
    ex. 6 = pass 2, 3
    ex. 1 = pass 0
*/
bool NtupleMaker::triggerDecision(edm::Handle<edm::TriggerResults> &hltR, int iTrigger) {
    return hltR->wasrun(iTrigger) && hltR->accept(iTrigger) && !hltR->error(iTrigger);
}

void NtupleMaker::analyzeTrigger(edm::Handle<edm::TriggerResults> &hltR, edm::Handle<trigger::TriggerEvent> &hltE, const string& triggerName) {

    const unsigned int n = hltConfig_.size();
    const unsigned int triggerIndex = hltConfig_.triggerIndex(triggerName);
    vout<<"n = "<<n<<" triggerIndex = "<<triggerIndex<<" size = "<<hltConfig_.size()<<endl;
    vout<<"Analyze triggerName: "<<triggerName<<endl;
    if(triggerIndex >= n) {
        vout<<"FourmuonAnalyzer4::analyzeTrigger: path "<<triggerName<<" - not found!"<<endl;
        return;
    }
    const unsigned int moduleIndex = hltR->index(triggerIndex);
    const unsigned int m = hltConfig_.size(triggerIndex);
    const vector<string>& moduleLabels = hltConfig_.moduleLabels(triggerIndex);
    if(Total_events_analyzed == 1){
        vout<<"FourmuonAnalyzer::analyzeTrigger: path "<<triggerName<<" ["<<triggerIndex<<"]"<<endl;
        vout<<"n = "<< n<<" triggerIndex = "<<triggerIndex<<" m = "<<m<<endl;
        vout<<"moduleLabels = "<<moduleLabels.size()<<" moduleIndex = "<<moduleIndex<<endl;
        // Results from TriggerResults product
        vout<<"Trigger path status:"<<" WasRun="<<hltR->wasrun(triggerIndex)<<" Accept="<<hltR->accept(triggerIndex)<<" Error ="<<hltR->error(triggerIndex)<<endl;
        vout<<"Last active module - label/type: "<<moduleLabels[moduleIndex]<<"/"<<hltConfig_.moduleType(moduleLabels[moduleIndex])<<" ["<<moduleIndex<<" out of 0-"<<(m - 1)<<" on this path]"<<endl;
    }
    assert(moduleIndex < m);
    // Results from TriggerEvent product - Attention: must look only for modules actually run in this path for this event!
    vector<GlobalVector> passMomenta;
    for(unsigned int j = 0; j <= moduleIndex; ++j) {
        const string& moduleLabel = moduleLabels[j];
        const string moduleType = hltConfig_.moduleType(moduleLabel);
        // check whether the module is packed up in TriggerEvent product
        const unsigned int filterIndex = hltE->filterIndex(edm::InputTag(moduleLabel, "", hltName_));
        if(Total_events_analyzed == 1) vout<<"j = "<<j<<" modLabel/moduleType = "<<moduleLabel<<"/"<<moduleType<<" filterIndex = "<<filterIndex<<" sizeF = "<<hltE->sizeFilters()<<endl;
        if(filterIndex < hltE->sizeFilters()) {
            if(Total_events_analyzed == 1) vout<<"'L3' (or 'L1', 'L2') filter in slot "<<j<<" - label/type "<<moduleLabel<<"/"<<moduleType<<endl;
            const Vids& VIDS(hltE->filterIds(filterIndex));
            const Keys& KEYS(hltE->filterKeys(filterIndex));
            const size_type nI(VIDS.size());
            const size_type nK(KEYS.size());
            assert(nI == nK);
            const size_type n(max(nI,nK));
            if(Total_events_analyzed == 1) vout<<n<<" accepted 'L3' (or 'L1', 'L2') objects found:"<<endl;
            const TriggerObjectCollection& TOC(hltE->getObjects());
            for (size_type i = 0; i != n; ++i) {
                if(i == 0) passMomenta.clear();
                const TriggerObject& TO(TOC[KEYS[i]]);
                GlobalVector momentumT0(TO.px(), TO.py(), TO.pz());
                if(Total_events_analyzed == 1) vout<<"i = "<<i<<" moduleLabel/moduleType : "<<moduleLabel<<"/"<<moduleType<<endl;
                if(Total_events_analyzed == 1) vout<<i<<" "<<VIDS[i]<<"/"<<KEYS[i]<<": "<<TO.id()<<" "<<TO.pt()<<" "<<TO.eta()<<" "<<TO.phi()<<" "<<TO.mass()<<endl;
                if(TO.id() != 13 && TO.id() != -13 && TO.id() != 0) continue;
                // TO.id() --> L1 Mu is always 0 (?)
                vout<<" current moduleType = "<<moduleType<<endl;
                // HLT_Dimuon0_Jpsi_Muon_v5
                if("HLTL1TSeed" == moduleType && "hltL1sTripleMu0orTripleMu500" == moduleLabel) vout<<" L1 object found"<<endl;
                // HLT_Dimuon0_Jpsi_Muon_v5
                else if("HLTMuonL2FromL1TPreFilter" == moduleType) vout<<" L2 object found"<<endl;
                else if("HLTMuonL3PreFilter" == moduleType || "HLTMuonIsoFilter" == moduleType) vout<<" L3 (highEff) object found"<<endl;
                else if("HLTDiMuonGlbTrkFilter" == moduleType) vout<<" L3 (lowEff) object found"<<endl;
                else if("HLT2MuonMuonDZ" == moduleType || "HLTMuonIsoFilter" == moduleType || ("HLTMuonL3PreFilter" == moduleType)) vout<<" HLT object found"<<endl;
                // HLT_Dimuon0_Jpsi_Muon_v5
                else if("HLTMuonDimuonL3Filter" == moduleType) vout<<" HLT L3 filter object found"<<endl;
                // HLT muons decayed from Upsilon candidate
                else if("HLTDisplacedmumuFilter" == moduleType) vout<<"HLT object from (Jpsi candidate) found"<<endl;
                else continue;
                passMomenta.push_back(momentumT0);
            }
        }
        if("HLTMuonL3PreFilter" == moduleType){
            for(unsigned int i = 0; i < passMomenta.size(); ++i) muonTriggeredVecs.push_back(passMomenta[i]);
            vout<<" HLT obj FOUND ; current size = "<<muonTriggeredVecs.size()<<endl;
        }
        if("HLTDisplacedmumuFilter" == moduleType){// HLT muons decayed from Upsilon candidate
            for(unsigned int i = 0; i < passMomenta.size(); ++i) JpsiTriggeredVecs.push_back(passMomenta[i]);
            vout<<"HLT object from (Upsilon candidate) found: current size ="<<JpsiTriggeredVecs.size()<<endl;
        }
        passMomenta.clear();
    }
    return;
}

bool NtupleMaker::findTrigger(edm::Handle<edm::TriggerResults> &hltR, vector<string>& triggersToCheck, vector<string>& triggerNameFound) {
    triggerNameFound.clear();
    vout<<"findTrigger()... "<<endl<<"Request: "<<endl;
    for(unsigned int iT = 0; iT < triggersToCheck.size(); ++iT) vout<<"name ["<<iT<<"] = "<<triggersToCheck[iT]<<endl;

    bool triggerFound = false;
    vector<string> hlNames = triggerNames.triggerNames();
    for(uint iT = 0; iT < hlNames.size(); ++iT) {
        if(Total_events_analyzed == 1) vout<<" iT = "<<iT<<" hlNames[iT] = "<<hlNames[iT]<<" : wasrun = "<<hltR->wasrun(iT)
        <<" accept = "<<hltR->accept(iT)<<" !error = "<<!hltR->error(iT)<<endl;
        for(uint imyT = 0; imyT < triggersToCheck.size(); ++imyT){
            if(string::npos != hlNames[iT].find(triggersToCheck[imyT])) {
                if(Total_events_analyzed == 1) vout<<" Trigger "<<hlNames[iT]<<" found to be compatible with the requested. "<<endl;
                triggerNameFound.push_back(hlNames[iT]);
                if(triggerDecision(hltR, iT)) triggerFound = true;
            }
        }
    }
    return triggerFound;
}

UInt_t NtupleMaker::getTriggerBits(const edm::Event& iEvent) {
	UInt_t itrigger = 0;
    // trigger collection
	edm::Handle<edm::TriggerResults> triggerResults_handle;
	iEvent.getByToken(triggerResultsTok_, triggerResults_handle);
    // if(triggerResults_handle.isValid()) { 
    //     string testTriggerName;
    //     const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
    //     for(unsigned int trigIndex = 0; trigIndex < TheTriggerNames.size(); trigIndex++) {
    //         testTriggerName = TheTriggerNames.triggerName(trigIndex);
    //         cout<<testTriggerName.c_str()<<endl;
    //     }
    // }
	if(triggerResults_handle.isValid()) {// trigger collection 
		const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
		vector <unsigned int> bits_0, bits_1, bits_2, bits_3, bits_4, bits_5, bits_6, bits_7, bits_8, bits_9;
		for (int version = 1; version<20; version++) {
			stringstream ss0,ss1,ss2,ss3,ss4,ss5,ss6,ss7,ss8,ss9;
			ss0<<"HLT_Dimuon25_Jpsi_v"<<version;
			//ss0<<"HLT_Dimuon16_Jpsi_v"<<version;
			bits_0.push_back(TheTriggerNames.triggerIndex(edm::InputTag(ss0.str()).label().c_str()));

			ss1<<"HLT_Dimuon0_Upsilon_Muon_v"<<version;
			//ss1<<"HLT_Dimuon13_PsiPrime_v"<<version;
			bits_1.push_back(TheTriggerNames.triggerIndex(edm::InputTag(ss1.str()).label().c_str()));

			ss2<<"HLT_Dimuon13_Upsilon_v"<<version;
			bits_2.push_back(TheTriggerNames.triggerIndex(edm::InputTag(ss2.str()).label().c_str()));

			ss3<<"HLT_Dimuon0_Jpsi_Muon_v"<<version; // 2016 charmonium
			//ss3<<"HLT_Dimuon10_Jpsi_Barrel_v"<<version;
			bits_3.push_back(TheTriggerNames.triggerIndex(edm::InputTag(ss3.str()).label().c_str()));

            ss4<<"HLT_Dimuon0_Jpsi3p5_Muon2_v"<<version; ///2017-18 charmonium
			//ss4<<"HLT_TripleMu5_v"<<version; 
			//ss4<<"HLT_Dimuon8_PsiPrime_Barrel_v"<<version;
			bits_4.push_back(TheTriggerNames.triggerIndex(edm::InputTag(ss4.str()).label().c_str()));

			ss5<<"HLT_Dimuon12_Upsilon_y1p4_v"<<version;
			//ss5<<"HLT_Dimuon8_Upsilon_Barrel_v"<<version;
			bits_5.push_back(TheTriggerNames.triggerIndex(edm::InputTag(ss5.str()).label().c_str()));

			ss6<<"HLT_Dimuon20_Jpsi_v"<<version;
			bits_6.push_back(TheTriggerNames.triggerIndex(edm::InputTag(ss6.str()).label().c_str()));

			ss7<<"HLT_Trimuon2_Upsilon5_Muon_v"<<version;
			//ss7<<"HLT_Dimuon14_Phi_Barrel_Seagulls_v"<<version;
			//ss7<<"HLT_Dimuon0_Phi_Barrel_v"<<version;
			bits_7.push_back(TheTriggerNames.triggerIndex(edm::InputTag(ss7.str()).label().c_str()));

			ss8<<"HLT_Trimuon5_3p5_2_Upsilon_Muon_v"<<version;
			//ss8<<"HLT_DoubleMu4_3_Bs_v"<<version;
			bits_8.push_back(TheTriggerNames.triggerIndex(edm::InputTag(ss8.str()).label().c_str()));

			ss9<<"HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_v"<<version;
			//ss9<<"HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v"<<version;
			bits_9.push_back(TheTriggerNames.triggerIndex(edm::InputTag(ss9.str()).label().c_str()));
		}
		for(unsigned int i=0; i<bits_0.size(); i++) {
			unsigned int bit = bits_0[i];
			if(bit < triggerResults_handle->size()){
				if(triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
					itrigger += 1;
					break;
				}
			}
		}
		for(unsigned int i=0; i<bits_1.size(); i++) {
			unsigned int bit = bits_1[i];
			if(bit < triggerResults_handle->size()){
				if(triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
					itrigger += 2;
					break;
				}
			}
		}
		for(unsigned int i=0; i<bits_2.size(); i++) {
			unsigned int bit = bits_2[i];
			if(bit < triggerResults_handle->size()){
				if(triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
					itrigger += 4;
					break;
				}
			}
		}
		for(unsigned int i=0; i<bits_3.size(); i++) {
			unsigned int bit = bits_3[i];
			if(bit < triggerResults_handle->size()){
				if(triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
					itrigger += 8;
					break;
				}
			}
		}
		for(unsigned int i=0; i<bits_4.size(); i++) {
			unsigned int bit = bits_4[i];
			if(bit < triggerResults_handle->size()){
				if(triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
					itrigger += 16;
					break;
				}
			}
		}
		for(unsigned int i=0; i<bits_5.size(); i++) {
			unsigned int bit = bits_5[i];
			if(bit < triggerResults_handle->size()){
				if(triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
					itrigger += 32;
					break;
				}
			}
		}
		for(unsigned int i=0; i<bits_6.size(); i++) {
			unsigned int bit = bits_6[i];
			if(bit < triggerResults_handle->size()){
				if(triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
					itrigger += 64;
					break;
				}
			}
		}
		for(unsigned int i=0; i<bits_7.size(); i++) {
			unsigned int bit = bits_7[i];
			if(bit < triggerResults_handle->size()){
				if(triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
					itrigger += 128;
					break;
				}
			}
		}
		for(unsigned int i=0; i<bits_8.size(); i++) {
			unsigned int bit = bits_8[i];
			if(bit < triggerResults_handle->size()){
				if(triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) { 
					itrigger += 256;
					break;
				}   
			}   
		}   
		for(unsigned int i=0; i<bits_9.size(); i++) {
			unsigned int bit = bits_9[i];
			if(bit < triggerResults_handle->size()){
				if(triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) { 
					itrigger += 512;
					break;
				}   
			}   
		} 
    } 
	return itrigger;
}

bool NtupleMaker::triggerMatch(TLorentzVector mu1p4, TLorentzVector mu2p4, TLorentzVector mu3p4, TLorentzVector mu4p4) {
    int nTrigMuon = 0;
    double drCut = 0.3, dpCut = 1, dpRatioCut = 0.1;
    // match J/psi muons
    vout<<"Trigger matching for dimuon candidate"<<endl;
    double reco1_eta = mu1p4.Eta(), reco1_phi = mu1p4.Phi(), reco1_pt = mu1p4.Pt();
    double reco2_eta = mu2p4.Eta(), reco2_phi = mu2p4.Phi(), reco2_pt = mu2p4.Pt();
    vout<<"This Dimuon candidate"<<" mu1pt:"<<reco1_pt<<" mu1eta:"<<reco1_eta<<" mu1phi:"<<reco1_phi<<endl;
    vout<<"This Dimuon candidate"<<" mu2pt:"<<reco2_pt<<" mu1eta:"<<reco2_eta<<" mu1phi:"<<reco2_phi<<endl;
    double dR1_minimum = drCut, dR2_minimum = drCut, dPt1 = dpCut, dPt2 = dpCut, dPt1OvP = dpRatioCut, dPt2OvP = dpRatioCut;
    double hlt_pt1, hlt_pt2;
    vout<<"JpsiTriggeredVecs.size();"<<JpsiTriggeredVecs.size()<<endl;
    for(uint iTrig = 0; iTrig < JpsiTriggeredVecs.size(); ++iTrig) {// loop over all HLT objects
        double hlt_pt = JpsiTriggeredVecs[iTrig].perp();
        double hlt_eta = JpsiTriggeredVecs[iTrig].eta();
        double hlt_phi = JpsiTriggeredVecs[iTrig].phi();
        double dR1 = deltaR(reco1_eta, reco1_phi, hlt_eta, hlt_phi);
        double dR2 = deltaR(reco2_eta, reco2_phi, hlt_eta, hlt_phi);
        vout<<"iTrig:"<<iTrig<<" iTrigpT:"<<hlt_pt<<" hlt_phi:"<<hlt_phi<<endl;
        vout<<"dR1:"<<dR1<<" dPt1:"<<dPt1<<endl;
        vout<<"dR2:"<<dR2<<" dPt2:"<<dPt2<<endl;
        if(dR1 < dR2) {// matching to muon1 
            // checking best matching object with muon 1 
            if(dR1 >= dR1_minimum) continue;
            dR1_minimum = dR1;
            dPt1 = abs(reco1_pt - hlt_pt);
            dPt1OvP = dPt1 / reco1_pt;
            if(dR1 < drCut && (dPt1 < dpCut || dPt1OvP < dpRatioCut)) {
                vout<<"Matching Vtx sucessfull muPt1 = "<<reco1_pt<<" trigPt = "<<hlt_pt<<" dR = "<<dR1<<endl;
                hlt_pt1 = hlt_pt;
                // temp_L1muPt.push_back(hlt_pt);
            } //Found matching to muon 1
        }else {// matching to muon2
            // checking best matching object with muon 2
            if(dR2 >= dR2_minimum) continue;
            dR2_minimum = dR2;
            dPt2 = abs(reco2_pt - hlt_pt);
            dPt2OvP = dPt2 / reco2_pt;
            if(dR2 < drCut && (dPt2 < dpCut || dPt2OvP < dpRatioCut)) {
                vout<<"Matching Vtx sucessfull muPt2 = "<<reco2_pt<<" trigPt = "<<hlt_pt<<" dR = "<<dR2<<endl;
                hlt_pt2 = hlt_pt;
                // temp_L1muPt.push_back(hlt_pt);
            }
        } 
    }
    if(dR1_minimum < drCut && (dPt1 < dpCut || dPt1OvP < dpRatioCut)) nTrigMuon++;
    if(dR2_minimum < drCut && (dPt2 < dpCut || dPt2OvP < dpRatioCut)) nTrigMuon++;
    vout<<" AllTrigMu = "<<nTrigMuon<<endl<<" good fourMu  run/lumi/event : "<<run<<"/"<<lumi<<"/"<<event<<endl;
    if(nTrigMuon >= 2) vout<<"Matching good "<<endl;
    if(nTrigMuon < 2) return false;
    // match non-J/psi muons
    vout<<"Trigger matching for Rest of muons candidates"<<endl;
    double reco3_eta = mu3p4.Eta(), reco3_phi = mu3p4.Phi(), reco3_pt = mu3p4.Pt();
    double reco4_eta = mu4p4.Eta(), reco4_phi = mu4p4.Phi(), reco4_pt = mu4p4.Pt();
    vout<<"First muon candidate"<<" mu3pt:"<<reco3_pt<<" mu3eta:"<<reco3_eta<<" mu3phi:"<<reco3_phi<<endl;
    vout<<"Second muon candidate"<<" mu4pt:"<<reco4_pt<<" mu4eta:"<<reco4_eta<<" mu4phi:"<<reco4_phi<<endl;
    double dR3_minimum = 99, dR4_minimum = 99, dPt3 = 999, dPt4 = 999, dPt3OvP = 999, dPt4OvP = 999;
    vout<<"muonTriggeredVecs.size():"<<muonTriggeredVecs.size()<<endl;
    for (uint iTrig = 0; iTrig < muonTriggeredVecs.size(); ++iTrig) {
        double hlt_pt = muonTriggeredVecs[iTrig].perp();
        if(hlt_pt == hlt_pt1 || hlt_pt == hlt_pt2) continue;
        double hlt_eta = muonTriggeredVecs[iTrig].eta();
        double hlt_phi = muonTriggeredVecs[iTrig].phi();
        double dR3 = deltaR(reco3_eta, reco3_phi, hlt_eta, hlt_phi);
        double dR4 = deltaR(reco4_eta, reco4_phi, hlt_eta, hlt_phi);
        vout<<"iTrig:"<<iTrig<<" iTrigpT:"<<hlt_pt<<" hlt_phi:"<<hlt_phi<<endl;
        vout<<"dR3:"<<dR3<<" dPt3:"<<dPt3<<endl;
        vout<<"dR4:"<<dR4<<" dPt4:"<<dPt4<<endl;
        if(dR3 < dR4) {
            if(dR3 >= dR3_minimum) continue;
            dR3_minimum = dR3;
            dPt3 = abs(reco3_pt - hlt_pt);
            dPt3OvP = dPt3 / reco3_pt;
            if(dR3 < drCut && (dPt3 < dpCut || dPt3OvP < dpRatioCut)) {
                vout<<"Matching L3 sucessfull muPt1 = "<<reco3_pt<<" trigPt = "<<hlt_pt<<" dR = "<<dR3<<endl;
                // temp_L1muPt.push_back(hlt_pt);
            }
        }else {
            if(dR4 >= dR4_minimum) continue;
            dR4_minimum = dR4;
            dPt4 = abs(reco4_pt - hlt_pt);
            dPt4OvP = dPt4 / reco4_pt;
            if(dR4 < drCut && (dPt4 < dpCut || dPt4OvP < dpRatioCut)) {
                vout<<"Matching L3 sucessfull muPt2 = "<<reco4_pt<<" trigPt = "<<hlt_pt<<" dR = "<<dR4<<endl;
                // temp_L1muPt.push_back(hlt_pt);
            }
        }
    }
    if(dR3_minimum < drCut && (dPt3 < dpCut || dPt3OvP < dpRatioCut)) nTrigMuon++;
    if(dR4_minimum < drCut && (dPt4 < dpCut || dPt4OvP < dpRatioCut)) nTrigMuon++;
    vout<<" All Rest Trigger matched muons size = "<<nTrigMuon<<endl<<" run/lumi/event : "<<run<<"/"<<lumi<<"/"<<event<<endl;
    if(nTrigMuon >= 3) vout<<"Atleast one of Rest muon Match with HLT Object "<<endl;
    return nTrigMuon >= 3;
}

// bool NtupleMaker::TriggerMatch_restMuons(TLorentzVector mu3p4, TLorentzVector mu4p4) {}
// ------------ method called for each event  ------------
void NtupleMaker::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace pat;
    
	edm::Handle<pat::CompositeCandidateCollection> dimuons;
	iEvent.getByToken(dimuon_Label, dimuons);

	// edm::Handle<pat::CompositeCandidateCollection> conversions;
	// iEvent.getByToken(conversion_Label, conversions);

	edm::Handle< edm::View<pat::Muon> > muons;
	iEvent.getByToken(muon_Label, muons);

	edm::Handle<reco::VertexCollection> primaryVertices_handle;
	iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);
	reco::Vertex thePrimaryV;

	edm::Handle<reco::BeamSpot> theBeamSpot;
	iEvent.getByToken(bs_Label,theBeamSpot);
	reco::BeamSpot bs = *theBeamSpot;

    reco::VertexCollection primaryVertices = *primaryVertices_handle;
    primaryVertices.push_back(reco::Vertex(bs.position(), bs.covariance3D()));
    REevt_x.clear();
    REevt_y.clear();
    REevt_z.clear();
    for(auto it = primaryVertices.begin(); it != primaryVertices.end(); ++it) {
        REevt_x.push_back(it->x());
        REevt_y.push_back(it->y());
        REevt_z.push_back(it->z());
    }
    thePrimaryV = reco::Vertex(*(primaryVertices.begin()));
	// if(primaryVertices_handle->begin() != primaryVertices_handle->end()) thePrimaryV = reco::Vertex(*(primaryVertices_handle->begin()));
	// else thePrimaryV = reco::Vertex(bs.position(), bs.covariance3D());
	/*
        edm::Handle<reco::ConversionCollection> conversionHandle;
        iEvent.getByLabel(conversion_Label,conversionHandle);

        edm::Handle<reco::PFCandidateCollection> pfcandidates;
        iEvent.getByLabel("particleFlow",pfcandidates);
        const reco::PFCandidateCollection pfphotons = selectPFPhotons(*pfcandidates);
	*/
    ++Total_events_analyzed;

    if(int(iEvent.id().run()) != runNumber){
        runNumber = iEvent.id().run();
        vout<<" New run : "<<iEvent.id().run()<<endl;
        const edm::Run * iRun_c = &iEvent.getRun();
        edm::Run * iRun = const_cast <edm::Run*> (iRun_c);
        beginRun(*iRun, iSetup);
    }

    edm::Handle<edm::TriggerResults> hltR;
    edm::Handle<trigger::TriggerEvent> hltE;
    iEvent.getByToken(triggerEventTok_,hltE);
    iEvent.getByToken(triggerResultsTok_, hltR);
    vector<string>  hlNames;
    triggerNames = iEvent.triggerNames(*hltR);
    hlNames=triggerNames.triggerNames();
    vector < string > triggersFoundToApply;
    bool theTriggerPassed = findTrigger(hltR, triggerList, triggersFoundToApply);
    if(theTriggerPassed){
        vout<<" theTriggerPassed = "<<theTriggerPassed<<" trigFound = "<<triggersFoundToApply.size()<<endl;
        ++Total_events_analyzed_triggered;
    }
    trigger = getTriggerBits(iEvent);
    run     = iEvent.id().run();
    lumi    = iEvent.id().luminosityBlock();
    event   = iEvent.id().event();
    vout<<"run: "<<run<<endl;
    vout<<"lumi: "<<lumi<<endl;
    vout<<"event: "<<event<<endl;

    // allL1TrigMuons.clear();
    // allL2TrigMuons.clear();
    // allL3TrigMuons.clear();
    // allMuL1TriggerVectors.clear();
    // allMuL2TriggerVectors.clear();
    // allMuL3TriggerVectors_lowEff.clear();
    // allMuL3TriggerVectors_highEff.clear();
    muonTriggeredVecs.clear();
    JpsiTriggeredVecs.clear();

    GEmu_pt.clear();
    GEmu_eta.clear();
    GEmu_phi.clear();
    GEmu_mass.clear();
    GEmu_recoId.clear();

    GEJpsi_recoId = -1;
    GEJpsi_matchGEN = false;
    GEJpsi_passID = false;
    GEJpsi_passVtx = false;
    GEpsi2S_recoId = -1;
    GEpsi2S_matchGEN = false;
    GEpsi2S_passID = false;
    GEpsi2S_passVtx = false;
    GEevt_passHLT = false;
    GEevt_passAll = false;
    GEevt_fourMuMass = 0.;

    REmu_pt.clear();
    REmu_eta.clear();
    REmu_phi.clear();
    REmu_mass.clear();
    REmu_isSoft.clear();
    REmu_passCut.clear();
    REJpsi_pt.clear();
    REJpsi_eta.clear();
    REJpsi_y.clear();
    REJpsi_phi.clear();
    REJpsi_mass.clear();
    REJpsi_massErr.clear();
    REJpsi_vtxProb.clear();
    REJpsi_ctau.clear();
    REJpsi_ctauErr.clear();
    REJpsi_LxyPV.clear();
    REJpsi_sigLxy.clear();
    REJpsi_muId1.clear();
    REJpsi_muId2.clear();
    REJpsi_passCut.clear();
    REpsi2S_pt.clear();
    REpsi2S_eta.clear();
    REpsi2S_y.clear();
    REpsi2S_phi.clear();
    REpsi2S_mass.clear();
    REpsi2S_massErr.clear();
    REpsi2S_vtxProb.clear();
    REpsi2S_cstrVtxProb.clear();
    REpsi2S_ctau.clear();
    REpsi2S_ctauErr.clear();
    REpsi2S_LxyPV.clear();
    REpsi2S_sigLxy.clear();
    REpsi2S_muId1.clear();
    REpsi2S_muId2.clear();
    REpsi2S_passCut.clear();
    REevt_fourMuMass.clear();
    REevt_massChisq.clear();
    REevt_vtxProb.clear();
    REevt_d.clear();
    REevt_JpsiId.clear();
    REevt_psi2SId.clear();
    REevt_passHLT.clear();
    REevt_matchTrg.clear();
    REevt_fourMuFit.clear();
    REevt_twoDimuFit.clear();

    vout<<"triggersFoundToApply.size()"<<triggersFoundToApply.size()<<endl;
    for(unsigned int iTrig = 0; iTrig < triggersFoundToApply.size(); ++iTrig) {
        vout<<"triggersFoundToApply.at(iTrig)"<<triggersFoundToApply.at(iTrig)<<endl;
        analyzeTrigger(hltR, hltE, triggersFoundToApply.at(iTrig));
    }
    vout<<"Trigger analyzed Finished"<<endl;

    edm::Handle<reco::GenParticleCollection> prunedgenParticles;
    iEvent.getByToken(prunedgenParticlesSrc_, prunedgenParticles);
    edm::Handle<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag>> primaryVtx;
    iEvent.getByToken(primaryVtxToken_, primaryVtx);
    GEevt_vtx.clear();
    GEevt_vtx.push_back(primaryVtx->X());
    GEevt_vtx.push_back(primaryVtx->Y());
    GEevt_vtx.push_back(primaryVtx->Z());
    //edm::Handle<pat::PackedGenParticle> packedgenParticles;
    //iEvent.getByToken(packedgenParticlesSrc_, packedgenParticles);

	if(prunedgenParticles.isValid()) {
        vout<<"setting gen variables"<<endl;    
        setGENVariables(prunedgenParticles, *primaryVtx);
	}

	edm::ESHandle<TransientTrackBuilder> theB; 
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 
	edm::ESHandle<MagneticField> bFieldHandle;
	iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
    RefCountedKinematicTree tree_ups_part1;

    fourMuonFit(muons, bFieldHandle, bs, thePrimaryV,primaryVertices_handle);
    
    onia_tree->Fill();
} //end analyze
// ------------ method called once each job just before starting event loop  ------------
void NtupleMaker::beginJob() {
    vout<<"The job Begins"<<endl;
}
// ------------ method called once each job just after ending the event loop  ------------
void NtupleMaker::endJob() {
	cout <<"Job End"<<endl;
	cout<<"Total_events_analyzed: "<<Total_events_analyzed<<endl;
	cout<<"Total_events_analyzed_triggered: "<<Total_events_analyzed_triggered<<endl;
	cout<<"Total_events_dimuon_trg_matched: "<<Total_events_dimuon_trg_matched<<endl;
    cout<<"Total_events_Fourmuon: " <<Total_events_Fourmuon<<endl;
    cout<<"Total_events_soft4mu: "<<nevent_soft4mu<<endl;
    cout<<"Total_events_nevent_4mu_pt: "<<nevent_4mu_pt<<endl;
    cout<<"Total_events_nevent_4mu_eta " <<nevent_4mu_eta<<endl;
}
// ------------ method called when ending the processing of a run  ------------
void NtupleMaker::beginRun(edm::Run & iRun, edm::EventSetup const& iSetup) {
    vout<<" New run..."<<endl;
    //--- m_l1GtUtils.getL1GtRunCache(run, iSetup, true, false);
    bool hltConfigChanged;
    bool test = hltConfig_.init(iRun, iSetup, hltName_, hltConfigChanged);
    if(test) {
        vout<<" hltConfig_.size() = "<<hltConfig_.size()<<endl;
        // check iftrigger name in (new) config
        if(triggerName_!="@") {// "@" means: analyze all triggers in config
            const unsigned int n(hltConfig_.size());
            const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
            vout<<" triggerIndex = "<<triggerIndex<<endl;
            if(triggerIndex >= n) {
                vout<<"HLTEventAnalyzerAOD::beginRun:"<<" TriggerName "<<triggerName_<<" not available in (new) config!"<<endl;
                vout<<"Available TriggerNames are: "<<endl;
                hltConfig_.dump("Triggers");
            }
        }else vout<<" Bad trigger name"<<endl;
    }else cout<<"beginRun:"<<" HLT config extraction failure with process name "<<hltName_<<endl;
}

// ------------ method called when starting to processes a run  ------------
void NtupleMaker::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void NtupleMaker::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void NtupleMaker::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void NtupleMaker::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	// The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even ifit is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}
// void NtupleMaker::setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles, edm::Handle<pat::PackedGenParticle> packedgenParticles){
void NtupleMaker::setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles, const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag>& pv) {
    vector<const reco::Candidate*> temp_Jpsi, temp_psi2S, temp_JpsiMom, temp_psi2SMom;
    vector<vector<const reco::Candidate*>> temp_JpsiDau, temp_psi2SDau;
    bool passAcc = true, valid = false;
    TLorentzVector temp_fourMuLV;
    for(auto p = prunedgenParticles->begin(); p != prunedgenParticles->end(); p++) {
        if((p->pdgId() != 443 && p->pdgId() != 100443) || p->status() != 2) continue;
        vector<const reco::Candidate*> temp_mu;
        for(int i = 0; i < (int)p->numberOfDaughters(); i++) {
            const reco::Candidate *p_dau = p->daughter(i);
            while(p_dau->status() != 1 && p_dau->daughter(0)) p_dau = p_dau->daughter(0);
            if(fabs(p_dau->pdgId()) == 13 && p_dau->status() == 1) temp_mu.push_back(p_dau);//if daughter of J/psi is muon
        }
        if((int)temp_mu.size() != 2) continue;
        if(p->pdgId() == 443) {// if the particle is J/psi
            temp_JpsiDau.push_back(temp_mu);
            const reco::Candidate *Jpsi_self = &*p;
            temp_Jpsi.push_back(Jpsi_self);
            while(Jpsi_self->mother(0)->pdgId() == 443) Jpsi_self = Jpsi_self->mother(0);
            temp_JpsiMom.push_back(Jpsi_self->mother(0));
        }else {// if the particle is psi(2S)
            temp_psi2SDau.push_back(temp_mu);
            const reco::Candidate *psi2S_self = &*p;
            temp_psi2S.push_back(psi2S_self);
            while(psi2S_self->mother(0)->pdgId() == 100443) psi2S_self = psi2S_self->mother(0);
            temp_psi2SMom.push_back(psi2S_self->mother(0));
        }
    }
    for(int i = 0; i < (int)temp_Jpsi.size(); i++) {
        for(int j = 0; j < (int)temp_psi2S.size(); j++) {
            // if(temp_JpsiMom[i] != temp_psi2SMom[j]) continue;// only works for SPS
            valid = true;
            for(int k = 0; k < 2; k++) {
                GEmu_pt.push_back(temp_JpsiDau[i][k]->pt());
                GEmu_eta.push_back(temp_JpsiDau[i][k]->eta());
                GEmu_phi.push_back(temp_JpsiDau[i][k]->phi());
                GEmu_mass.push_back(temp_JpsiDau[i][k]->mass());
            }
            for(int k = 0; k < 2; k++) {
                GEmu_pt.push_back(temp_psi2SDau[j][k]->pt());
                GEmu_eta.push_back(temp_psi2SDau[j][k]->eta());
                GEmu_phi.push_back(temp_psi2SDau[j][k]->phi());
                GEmu_mass.push_back(temp_psi2SDau[j][k]->mass());
            }
            GEJpsi_pt = temp_Jpsi[i]->pt();
            GEJpsi_eta = temp_Jpsi[i]->eta();
            GEJpsi_y = temp_Jpsi[i]->rapidity();
            GEJpsi_phi = temp_Jpsi[i]->phi();
            GEJpsi_mass = temp_Jpsi[i]->mass();
            GEJpsi_ctau = calCtau(*(temp_Jpsi[i]), pv);
            GEJpsi_muId1 = 0;
            GEJpsi_muId2 = 1;
            GEpsi2S_pt = temp_psi2S[j]->pt();
            GEpsi2S_eta = temp_psi2S[j]->eta();
            GEpsi2S_y = temp_psi2S[j]->rapidity();
            GEpsi2S_phi = temp_psi2S[j]->phi();
            GEpsi2S_mass = temp_psi2S[j]->mass();
            GEpsi2S_ctau = calCtau(*(temp_psi2S[i]), pv);
            GEpsi2S_muId1 = 2;
            GEpsi2S_muId2 = 3;
            break;
        }
        if(valid) break;
    }
    for(int i = 0; i < (int)GEmu_pt.size(); i++) {
        GEmu_recoId.push_back(-1);// dummy
        TLorentzVector temp_muLV;
        temp_muLV.SetPtEtaPhiM(GEmu_pt[i], GEmu_eta[i], GEmu_phi[i], GEmu_mass[i]);
        temp_fourMuLV += temp_muLV;
        if(fabs(GEmu_eta[i]) > 2.4 || GEmu_pt[i] < 3.5) passAcc = false;
    }
    GEevt_valid = valid;
    if(GEevt_valid) GEevt_fourMuMass = temp_fourMuLV.M();
    GEevt_passAcc = passAcc;
    vout<<"Find "<<GEmu_pt.size()<<" muons in GEN level."<<endl;
    return;
}// finished Gen variables

void NtupleMaker::fourMuonFit(
    edm::Handle<edm::View<pat::Muon>> muons,
    edm::ESHandle<MagneticField> bFieldHandle,
    reco::BeamSpot bs,
    reco::Vertex thePrimaryV,
    edm::Handle<reco::VertexCollection> primaryVertices_handle) {
    // muon pre-selection
    vector<pat::Muon> selMuons;
    for(edm::View<pat::Muon>::const_iterator recoMuon = muons->begin(); recoMuon != muons->end(); ++recoMuon) {
        if(!(recoMuon->isGlobalMuon() || recoMuon->isTrackerMuon())) continue;
        if(recoMuon->pt() < 2.0 || fabs(recoMuon->eta()) > 2.4)  continue;
        selMuons.push_back(*recoMuon);
        REmu_pt.push_back(recoMuon->pt());
        REmu_eta.push_back(recoMuon->eta());
        REmu_phi.push_back(recoMuon->phi());
        REmu_mass.push_back(recoMuon->mass());
        bool isSoft = muon::isSoftMuon(*recoMuon, thePrimaryV);
        REmu_isSoft.push_back(isSoft);
        REmu_passCut.push_back(fabs(recoMuon->eta()) < 2.4 && recoMuon->pt() > 3.6 && isSoft);
    }
    int nSelMuon = selMuons.size();
    vout<<"Select ["<<nSelMuon<<"] global or tracker muons."<<endl;
    // GEN matching of muons
    vector<vector<double>> drMap;// drMap[GENid][RECOid]
    for(int i = 0; i < (int)GEmu_pt.size(); i++) {
        vector<double> drMapRow;
        for(int j = 0; j < nSelMuon; j++) drMapRow.push_back(deltaR(REmu_eta[j], REmu_phi[j], GEmu_eta[i], GEmu_phi[i]));
        drMap.push_back(drMapRow);
    }
    for(int i = 0; i < (int)GEmu_pt.size(); i++) {
        for(int j = 0; j < nSelMuon; j++) vout<<drMap[i][j]<<' ';
	vout<<endl;
    }
    int nMatchMuon = 0;
    for(int i = 0; i < (int)GEmu_pt.size(); i++) {
        double minDr = 0.03;
        int genMuonId = -1, recoMuonId = -1;
        for(int j = 0; j < (int)drMap.size(); j++) {
            for(int k = 0; k < nSelMuon; k++) {
                if(drMap[j][k] < minDr) {
                    minDr = drMap[j][k];
                    genMuonId = j;
                    recoMuonId = k;
                }
            }
        }
        if(minDr == 0.03) break;
        GEmu_recoId[genMuonId] = recoMuonId;
        // guarantee paired GEN/RECO muon won't be chosen again.
        for(int j = 0; j < (int)drMap.size(); j++) drMap[j][recoMuonId] = 1;
        for(int k = 0; k < nSelMuon; k++) drMap[genMuonId][k] = 1;
        nMatchMuon++;
    }
    vout<<"Find ["<<nMatchMuon<<"] Gen Matched muons."<<endl;
    GEJpsi_matchGEN = GEevt_valid && GEmu_recoId[GEJpsi_muId1] != -1 && GEmu_recoId[GEJpsi_muId2] != -1;
    GEpsi2S_matchGEN = GEevt_valid && GEmu_recoId[GEpsi2S_muId1] != -1 && GEmu_recoId[GEpsi2S_muId2] != -1;
    // single muon cut
    GEJpsi_passID = GEJpsi_matchGEN && REmu_passCut[GEmu_recoId[GEJpsi_muId1]] && REmu_passCut[GEmu_recoId[GEJpsi_muId2]];
    GEpsi2S_passID = GEpsi2S_matchGEN && REmu_passCut[GEmu_recoId[GEpsi2S_muId1]] && REmu_passCut[GEmu_recoId[GEpsi2S_muId2]];
    // combinition of muons
    KinematicConstraint *JpsiMassCons = new MassKinematicConstraint(JpsiMass, JpsiSigma), *psi2SMassCons = new MassKinematicConstraint(psi2SMass, psi2SSigma);
    vector<RefCountedKinematicTree> JpsiVtxTrees, psi2SVtxTrees;// Save J/psi,psi(2S) vertex tree
    for(int i = 0; i < nSelMuon; i++) {
        if(!REmu_passCut[i]) continue;
        for(int j = i + 1; j < nSelMuon; j++) {
            if(!REmu_passCut[j]) continue;
            if(selMuons[i].charge() + selMuons[j].charge()) continue;
            // select candidates for vertex fit
            TLorentzVector selMuoniLV, selMuonjLV;
            selMuoniLV.SetPtEtaPhiM(REmu_pt[i], REmu_eta[i], REmu_phi[i], REmu_mass[i]);
            selMuonjLV.SetPtEtaPhiM(REmu_pt[j], REmu_eta[j], REmu_phi[j], REmu_mass[j]);
            TLorentzVector selDimuonLV = selMuoniLV + selMuonjLV;
            if(selDimuonLV.M() < 2.7 || selDimuonLV.M() > 4.1) continue;
            // vertex fit without mass constraint
            reco::TrackRef muonTracki = selMuons[i].track(), muonTrackj = selMuons[j].track();
            reco::TransientTrack muonTrsTracki(muonTracki, &(*bFieldHandle)), muonTrsTrackj(muonTrackj, &(*bFieldHandle));
            KinematicParticleFactoryFromTransientTrack factory;
            vector<RefCountedKinematicParticle> fitMuons;
            fitMuons.push_back(factory.particle(muonTrsTracki, muonMass, float(0), float(0), muonSigma));
            fitMuons.push_back(factory.particle(muonTrsTrackj, muonMass, float(0), float(0), muonSigma));
            KinematicParticleVertexFitter dimuonFitter;
            RefCountedKinematicTree dimuonVtxFitTree = dimuonFitter.fit(fitMuons);
            vout<<"Fit process finish\n";
            if(dimuonVtxFitTree->isEmpty()) continue;
            vout<<"This dimuon candidate has valid vertex fit"<<endl;
            dimuonVtxFitTree->movePointerToTheTop();
            RefCountedKinematicParticle vFitDimu = dimuonVtxFitTree->currentParticle();
            KinematicState vFitDimuState = vFitDimu->currentState();
            RefCountedKinematicVertex vFitVertex = dimuonVtxFitTree->currentDecayVertex();
            Double_t dimuonMass = vFitDimuState.mass();
            Double_t dimuonMassErr = sqrt(vFitDimuState.kinematicParametersError().matrix()(6,6));   
            Double_t dimuonVtxProb = ChiSquaredProbability((double)(vFitVertex->chiSquared()), (double)(vFitVertex->degreesOfFreedom()));
            vout<<"Vertex fit probability: "<<dimuonVtxProb<<endl;
            TLorentzVector temp_dimuLV;
            temp_dimuLV.SetPxPyPzE(vFitDimuState.globalMomentum().x(), vFitDimuState.globalMomentum().y(), vFitDimuState.globalMomentum().z(), vFitDimuState.kinematicParameters().energy());
            Double_t ctau, ctauErr, LxyPV, sigLxy;
            calCtau(vFitVertex, vFitDimu, thePrimaryV, ctau, ctauErr, LxyPV, sigLxy);
            int muId1 = i, muId2 = i;
            if(REmu_pt[i] > REmu_pt[j]) muId2 = j;
            else muId1 = j;
            bool dimuCutPartial = dimuonVtxProb > 0.005 && temp_dimuLV.Pt() > 10 && fabs(temp_dimuLV.Rapidity()) < 2;
            // dimuon vertex fit with mass constraint
            if(dimuonMass > 2.7 && dimuonMass < 3.5) {
                KinematicParticleFitter massConstraintFitter;
                RefCountedKinematicTree JpsiVtxFitTree = massConstraintFitter.fit(JpsiMassCons, dimuonVtxFitTree);
                Double_t JpsiCstrVtxProb = 0;
                if(JpsiVtxFitTree->isValid()) {
                    JpsiVtxFitTree->movePointerToTheTop();
                    RefCountedKinematicVertex vFitCstrJpsiVtx = JpsiVtxFitTree->currentDecayVertex();
                    JpsiCstrVtxProb = ChiSquaredProbability((double)(vFitCstrJpsiVtx->chiSquared()), (double)(vFitCstrJpsiVtx->degreesOfFreedom()));
                    vout<<"J/psi mass constrained vertex fit probability: "<<JpsiCstrVtxProb<<endl;
                }
                REJpsi_pt.push_back(temp_dimuLV.Pt());
                REJpsi_eta.push_back(temp_dimuLV.Eta());
                REJpsi_y.push_back(temp_dimuLV.Rapidity());
                REJpsi_phi.push_back(temp_dimuLV.Phi());
                REJpsi_mass.push_back(dimuonMass);
                REJpsi_massErr.push_back(dimuonMassErr);
                REJpsi_vtxProb.push_back(dimuonVtxProb);
                REJpsi_cstrVtxProb.push_back(JpsiCstrVtxProb);
                REJpsi_ctau.push_back(ctau);
                REJpsi_ctauErr.push_back(ctauErr);
                REJpsi_LxyPV.push_back(LxyPV);
                REJpsi_sigLxy.push_back(sigLxy);
                REJpsi_muId1.push_back(muId1);
                REJpsi_muId2.push_back(muId2);
                REJpsi_passCut.push_back(dimuCutPartial && dimuonMass > 2.95 && dimuonMass < 3.25);
                if(GEJpsi_passID && ((GEmu_recoId[GEJpsi_muId1] == i && GEmu_recoId[GEJpsi_muId2] == j) || (GEmu_recoId[GEJpsi_muId1] == j && GEmu_recoId[GEJpsi_muId2] == i))) {
                    GEJpsi_recoId = REJpsi_pt.size() - 1;
                    GEJpsi_passVtx = REJpsi_passCut[GEJpsi_recoId];
                }
                JpsiVtxTrees.push_back(dimuonVtxFitTree);// Save J/psi vertex tree
            }
            if(dimuonMass > 3.3 && dimuonMass < 4.1) {
                KinematicParticleFitter massConstraintFitter;
                RefCountedKinematicTree psi2SVtxFitTree = massConstraintFitter.fit(psi2SMassCons, dimuonVtxFitTree);
                Double_t psi2SCstrVtxProb = 0;
                if(psi2SVtxFitTree->isValid()) {
                    psi2SVtxFitTree->movePointerToTheTop();
                    RefCountedKinematicVertex vFitCstrpsi2SVtx = psi2SVtxFitTree->currentDecayVertex();
                    psi2SCstrVtxProb = ChiSquaredProbability((double)(vFitCstrpsi2SVtx->chiSquared()), (double)(vFitCstrpsi2SVtx->degreesOfFreedom()));
                    vout<<"psi(2S) mass constrained vertex fit probability: "<<psi2SCstrVtxProb<<endl;
                }
                REpsi2S_pt.push_back(temp_dimuLV.Pt());
                REpsi2S_eta.push_back(temp_dimuLV.Eta());
                REpsi2S_y.push_back(temp_dimuLV.Rapidity());
                REpsi2S_phi.push_back(temp_dimuLV.Phi());
                REpsi2S_mass.push_back(dimuonMass);
                REpsi2S_massErr.push_back(dimuonMassErr);
                REpsi2S_vtxProb.push_back(dimuonVtxProb);
                REpsi2S_cstrVtxProb.push_back(psi2SCstrVtxProb);
                REpsi2S_ctau.push_back(ctau);
                REpsi2S_ctauErr.push_back(ctauErr);
                REpsi2S_LxyPV.push_back(LxyPV);
                REpsi2S_sigLxy.push_back(sigLxy);
                REpsi2S_muId1.push_back(muId1);
                REpsi2S_muId2.push_back(muId2);
                REpsi2S_passCut.push_back(dimuCutPartial && dimuonMass > 3.35 && dimuonMass < 4.05);
                if(GEpsi2S_passID && ((GEmu_recoId[GEpsi2S_muId1] == i && GEmu_recoId[GEpsi2S_muId2] == j) || (GEmu_recoId[GEpsi2S_muId1] == j && GEmu_recoId[GEpsi2S_muId2] == i))) {
                    GEpsi2S_recoId = REpsi2S_pt.size() - 1;
                    GEpsi2S_passVtx = REpsi2S_passCut[GEpsi2S_recoId];
                }
                psi2SVtxTrees.push_back(dimuonVtxFitTree);// Save psi(2S) vertex tree
            }
        }
    }
    int nJpsi = REJpsi_pt.size(), npsi2S = REpsi2S_pt.size();
    vout<<"In this event: ["<<nJpsi<<"] J/psi, ["<<npsi2S<<"] psi(2S)"<<endl;
    if(!nJpsi || !npsi2S) return;
    // check HLT
    bool evtPassHLT = trigger & 16;
    if (evtPassHLT) vout<<"Event passed trigger selection"<<endl;
    else return;
    GEevt_passHLT = GEJpsi_passVtx && GEpsi2S_passVtx && evtPassHLT;
    // four muon vertex fit and trigger matching
    for(int i = 0; i < nJpsi; i++) {
        if(!REJpsi_passCut[i]) continue;
        int i1 = REJpsi_muId1[i], i2 = REJpsi_muId2[i];
        for(int j = 0; j < npsi2S; j++) {
            if(!REpsi2S_passCut[j]) continue;
            int j1 = REpsi2S_muId1[j], j2 = REpsi2S_muId2[j];
            if(i1 == j1 || i1 == j2 || i2 == j1 || i2 == j2) continue;
            RefCountedKinematicTree JpsiVtxTree = JpsiVtxTrees[i], psi2SVtxTree = psi2SVtxTrees[j];
            JpsiVtxTree->movePointerToTheTop();
            psi2SVtxTree->movePointerToTheTop();
            // trigger matching
            TLorentzVector muoni1LV, muoni2LV, muonj1LV, muonj2LV;
            muoni1LV.SetPtEtaPhiM(REmu_pt[i1], REmu_eta[i1], REmu_phi[i1], REmu_mass[i1]);
            muoni2LV.SetPtEtaPhiM(REmu_pt[i2], REmu_eta[i2], REmu_phi[i2], REmu_mass[i2]);
            muonj1LV.SetPtEtaPhiM(REmu_pt[j1], REmu_eta[j1], REmu_phi[j1], REmu_mass[j1]);
            muonj2LV.SetPtEtaPhiM(REmu_pt[j2], REmu_eta[j2], REmu_phi[j2], REmu_mass[j2]);
            bool evtTrgMatch = triggerMatch(muoni1LV, muoni2LV, muonj1LV, muonj2LV);
            // Four muon vertex fit, 2 methods provided
            // Method 1: fit 4 muons altogether
            vector<RefCountedKinematicParticle> fitMuons;
            reco::TransientTrack muoni1TrsTrack(selMuons[i1].track(), &(*bFieldHandle));
            reco::TransientTrack muoni2TrsTrack(selMuons[i2].track(), &(*bFieldHandle));
            reco::TransientTrack muonj1TrsTrack(selMuons[j1].track(), &(*bFieldHandle));
            reco::TransientTrack muonj2TrsTrack(selMuons[j2].track(), &(*bFieldHandle));
            KinematicParticleFactoryFromTransientTrack factory;
            fitMuons.push_back(factory.particle(muoni1TrsTrack, muonMass, 0., 0., muonSigma));
            fitMuons.push_back(factory.particle(muoni2TrsTrack, muonMass, 0., 0., muonSigma));
            fitMuons.push_back(factory.particle(muonj1TrsTrack, muonMass, 0., 0., muonSigma));
            fitMuons.push_back(factory.particle(muonj2TrsTrack, muonMass, 0., 0., muonSigma));
            KinematicParticleVertexFitter fourMuFitter;
            RefCountedKinematicTree fourMuVtxFitTree = fourMuFitter.fit(fitMuons);
            Double_t fourMuMass = 0, fourMuVtxProb = -1;
            bool passFitFourMu = !fourMuVtxFitTree->isEmpty();
            if(passFitFourMu) {
                fourMuVtxFitTree->movePointerToTheTop();
                RefCountedKinematicParticle fourMu = fourMuVtxFitTree->currentParticle();
                if(fourMu->currentState().isValid()) fourMuMass = fourMu->currentState().mass();
                RefCountedKinematicVertex fourMuVtx = fourMuVtxFitTree->currentDecayVertex();
                if(fourMuVtx->vertexIsValid()) fourMuVtxProb = ChiSquaredProbability((double)(fourMuVtx->chiSquared()), (double)(fourMuVtx->degreesOfFreedom()));
            }
            vout<<"Four-muon vertex probability: "<<fourMuVtxProb<<endl;
            // // Method 2: split into 2 dimuons
            vector<RefCountedKinematicParticle> fitDimus;
            RefCountedKinematicParticle JpsiVtx = JpsiVtxTree->currentParticle(), psi2SVtx = psi2SVtxTree->currentParticle();
            fitDimus.push_back(JpsiVtx);
            fitDimus.push_back(psi2SVtx);
            KinematicParticleVertexFitter twoDimuFitter;
            RefCountedKinematicTree twoDimuVtxFitTree = twoDimuFitter.fit(fitDimus);
            bool passFitDimu = !twoDimuVtxFitTree->isEmpty();
            // Handel output variables
            Double_t d = calD(JpsiVtxTree->currentDecayVertex(), psi2SVtxTree->currentDecayVertex());
            Double_t massChisq = SQUARE((REJpsi_mass[i] - JpsiMass) / REJpsi_massErr[i]) + SQUARE((REpsi2S_mass[j] - psi2SMass) / REpsi2S_massErr[j]);
            int pos = lower_bound(REevt_massChisq.begin(), REevt_massChisq.end(), massChisq) - REevt_massChisq.begin();
            REevt_fourMuMass.insert(REevt_fourMuMass.begin() + pos, fourMuMass);
            REevt_massChisq.insert(REevt_massChisq.begin() + pos, massChisq);
            REevt_vtxProb.insert(REevt_vtxProb.begin() + pos, fourMuVtxProb);
            REevt_d.insert(REevt_d.begin() + pos, d);
            REevt_JpsiId.insert(REevt_JpsiId.begin() + pos, i);
            REevt_psi2SId.insert(REevt_psi2SId.begin() + pos, j);
            REevt_passHLT.insert(REevt_passHLT.begin() + pos, evtPassHLT);
            REevt_matchTrg.insert(REevt_matchTrg.begin() + pos, evtTrgMatch);
            REevt_fourMuFit.insert(REevt_fourMuFit.begin() + pos, passFitFourMu);
            REevt_twoDimuFit.insert(REevt_twoDimuFit.begin() + pos, passFitDimu);
            if(GEJpsi_recoId == i && GEpsi2S_recoId == j) GEevt_passAll = GEevt_passHLT && evtTrgMatch && passFitFourMu && (fourMuMass > 7.5);
        }
    }
    return;
}

void NtupleMaker::calCtau(RefCountedKinematicVertex& decayVrtx, RefCountedKinematicParticle& kinePart, Vertex& bs,
    Double_t& ctau, Double_t& ctauErr, Double_t& LxyPV, Double_t& sigLxy) {
    TVector3 vtx(decayVrtx->position().x(), decayVrtx->position().y(), 0);
    TVector3 pvtx(bs.position().x(), bs.position().y(), 0);
    TVector3 pperp(kinePart->currentState().globalMomentum().x(), kinePart->currentState().globalMomentum().y(), 0);
    TVector3 vdiff = vtx - pvtx;
    LxyPV = vdiff.Dot(pperp) / pperp.Mag();
    ctau = LxyPV * kinePart->currentState().mass() / pperp.Perp();

    GlobalError DecayErr = decayVrtx->error();
    GlobalError PrimaryErr = bs.error();
    AlgebraicVector3 vdiffp;
    vdiffp[0] = vdiff.x();
    vdiffp[1] = vdiff.y();
    vdiffp[2] = 0;
    double LxyErr2 = ROOT::Math::Similarity((AlgebraicSymMatrix33)(DecayErr.matrix() + PrimaryErr.matrix()), vdiffp);
    double cosAlpha = vdiff.Dot(pperp) / (vdiff.Perp() * pperp.Perp());
    sigLxy = LxyPV / (TMath::Sqrt(LxyErr2) / vdiff.Perp() * cosAlpha);
    AlgebraicVector vpperp(3);
    vpperp[0] = pperp.x();
    vpperp[1] = pperp.y();
    vpperp[2] = 0.;
    AlgebraicSymMatrix vXYe = asHepMatrix(DecayErr.matrix()) + asHepMatrix(PrimaryErr.matrix());
    ctauErr = sqrt(vXYe.similarity(vpperp)) * kinePart->currentState().mass() / (pperp.Perp2());
    
    return;
}
double NtupleMaker::calCtau(const reco::Candidate& p, const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag>& pV) {
    TVector3 vtx(p.vx(), p.vy(), 0);
    // cout<<"sec: ("<<p.vx()<<", "<<p.vy()<<")\n";
    TVector3 pvtx(pV.X(), pV.Y(), 0);
    // cout<<"pri: ("<<pV.X()<<", "<<pV.Y()<<")\n";
    TVector3 pperp(p.px(), p.py(), 0);
    TVector3 vdiff = vtx - pvtx;
    double LxyPV = vdiff.Dot(pperp) / pperp.Mag();
    return LxyPV * p.mass() / pperp.Perp();
}

Double_t NtupleMaker::calD(RefCountedKinematicVertex decayVrtx1, RefCountedKinematicVertex decayVrtx2) {
    TVector3 vperp;
    vperp.SetXYZ(decayVrtx1->position().x() - decayVrtx2->position().x(), decayVrtx1->position().y() - decayVrtx2->position().y(), 0);
    double dist = vperp.Perp();
    GlobalError DecayErr1 = decayVrtx1->error();
    GlobalError DecayErr2 = decayVrtx2->error();
    AlgebraicVector3 aperp;
    aperp[0] = vperp.x();
    aperp[1] = vperp.y();
    aperp[2] = 0;
    double DistanceErr2 = ROOT::Math::Similarity((AlgebraicSymMatrix33)(DecayErr1.matrix() + DecayErr2.matrix()), aperp);
    return dist / (TMath::Sqrt(DistanceErr2) / dist);
}

// define this as a plug-in
DEFINE_FWK_MODULE(NtupleMaker);
