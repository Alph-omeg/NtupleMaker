// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      NtupleMaker
// 
// Description: Analyze pp -> J/psi+psi(2S) -> mu+mu-mu+mu-
//
// Author:  Zhen Hu
//
/*
    Version: Run II 2018 MINIAOD
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

// #define MONTE_CARLO// obselete
// #define VERBOSE// switch between verbose output and concise output
#ifdef VERBOSE
    #define vout cout
#else
    #define vout if(false)cout
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
const float JpsiSigma = 0.00004;
enum TrgLabel {L1, L3, Vtx};
//
// class declaration
//
class NtupleMaker: public edm::EDAnalyzer {
    public:
        explicit NtupleMaker(const edm::ParameterSet &);
        ~NtupleMaker();
        static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

    private:
        bool triggerMatch(const edm::Event &iEvent, const edm::Handle<edm::TriggerResults>& hltR, pat::Muon& muon, enum TrgLabel flag);
        void setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles);
        void fourMuonFit(const edm::Event &iEvent, const edm::Handle<edm::TriggerResults>& hltR, edm::Handle<edm::View<pat::Muon>> muons, edm::Handle<edm::View<pat::PackedCandidate>> candidates, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV, edm::Handle<reco::VertexCollection>);
        void calCtau(RefCountedKinematicVertex& decayVrtx, RefCountedKinematicParticle& kinePart, Vertex& bs, Double_t& ctau, Double_t& ctauErr, Double_t& LxyPV, Double_t& sigLxy);
        Double_t calD(RefCountedKinematicVertex decayVrtx1, RefCountedKinematicVertex decayVrtx2);

        virtual void beginJob();
        virtual void analyze(const edm::Event &, const edm::EventSetup &);
        virtual void endJob();

        virtual void beginRun(edm::Run &, edm::EventSetup const&);
        virtual void endRun(edm::Run const &, edm::EventSetup const &);
        virtual void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);
        virtual void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);
		// ----------member data ---------------------------
        // edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_Label;
		// edm::EDGetTokenT<pat::CompositeCandidateCollection> conversion_Label;
		edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
		edm::EDGetTokenT<reco::BeamSpot> bs_Label;
		edm::EDGetTokenT<edm::View<pat::Muon>> muon_Label;
		edm::EDGetTokenT<edm::TriggerResults> triggerResultsTok_;
        edm::EDGetTokenT<trigger::TriggerEvent> triggerEventTok_;
        edm::EDGetTokenT<edm::View<pat::PackedCandidate>> candidate_Label;
        //edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
        
        // GEN level muon
        vector<Double_t> GEmu_pt, GEmu_eta, GEmu_phi, GEmu_mass;
        vector<int> GEmu_recoId;
        // GEN level J/psi_1
        Double_t GEJpsi1_pt, GEJpsi1_eta, GEJpsi1_y, GEJpsi1_phi, GEJpsi1_mass;
        int GEJpsi1_recoId, GEJpsi1_muId1, GEJpsi1_muId2;
        bool GEJpsi1_matchGEN, GEJpsi1_passID, GEJpsi1_passVtx;
        // GEN level J/psi_2
        Double_t GEJpsi2_pt, GEJpsi2_eta, GEJpsi2_y, GEJpsi2_phi, GEJpsi2_mass;
        int GEJpsi2_recoId, GEJpsi2_muId1, GEJpsi2_muId2;
        bool GEJpsi2_matchGEN, GEJpsi2_passID, GEJpsi2_passVtx;
        // GEN level event
        Double_t GEevt_fourMuMass, GEevt_muPtMax;
        bool GEevt_valid, GEevt_passAcc, GEevt_passHLT, GEevt_matchTrg;
        // RECO level muon
        vector<Double_t> REmu_pt, REmu_eta, REmu_phi, REmu_mass;
        vector<bool> REmu_isSoft, REmu_passCut;
        vector<int> REmu_pvAsc;

        // test functionality: muon candidates
        vector<Double_t> TEmu_pt, TEmu_eta, TEmu_phi, TEmu_mass;
        vector<int> TEmu_pvAsc;

        // RECO level J/psi
        vector<Double_t> REJpsi_pt, REJpsi_eta, REJpsi_y, REJpsi_phi, REJpsi_mass, REJpsi_massErr, REJpsi_vtxProb, REJpsi_cstrVtxProb,
        REJpsi_ctau, REJpsi_ctauErr, REJpsi_LxyPV, REJpsi_sigLxy, REJpsi_X, REJpsi_Y, REJpsi_Z;
        vector<int> REJpsi_muId1, REJpsi_muId2;
        vector<bool> REJpsi_passCut;
        // RECO level event
        vector<double> REvtx_X, REvtx_Y, REvtx_Z;
        vector<Double_t> REevt_fourMuMass, REevt_massChisq, REevt_vtxProb, REevt_L1muPtMax, REevt_d, REevt_X, REevt_Y, REevt_Z;
        vector<int> REevt_JpsiId1, REevt_JpsiId2;
        vector<bool> REevt_passHLT, REevt_matchTrg, REevt_fourMuFit, REevt_twoDimuFit, REevt_samePV;

        bool evtPassHLT;
        int Total_events_analyzed;
        int Total_events_triggered;
        int Total_events_trg_matched;
        // some counters for cut flow
		ULong64_t run;
		ULong64_t lumi;
        ULong64_t event;
 
        string hltName_;
        string triggerName_;
        
        vector<string> triggerList;
        HLTConfigProvider hltConfig_;
		TTree *onia_tree;

        edm::EDGetTokenT<reco::GenParticleCollection> prunedgenParticlesSrc_;
};
//
// constructors and destructor
//
NtupleMaker::NtupleMaker(const edm::ParameterSet & iConfig):
    // dimuon_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("dimuons"))),
    // conversion_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("conversions"))),
    primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
    bs_Label(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("offlineBeamSpot"))),
    muon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter< edm::InputTag>("muons"))),
    triggerResultsTok_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
    // triggerEventTok_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("TriggerSummaryAOD"))),
    candidate_Label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Candidates"))),
    triggerList(iConfig.getUntrackedParameter<vector<string>>("triggerList")) {
	
    edm::Service <TFileService> fs;
	onia_tree = fs->make<TTree>("oniaTree", "Tree of FourMuon");

    Total_events_analyzed = 0;
    Total_events_triggered = 0 ;
    Total_events_trg_matched = 0;

    hltName_ = "HLT";
    triggerName_ = "@"; // "@" means: analyze all triggers in config

    onia_tree->Branch("run", &run, "run/l");
    onia_tree->Branch("lumi", &lumi, "lumi/l");
    onia_tree->Branch("event", &event, "event/l");

    // GEN level muon
    onia_tree->Branch("GEmu_pt", &GEmu_pt);
    onia_tree->Branch("GEmu_eta", &GEmu_eta);
    onia_tree->Branch("GEmu_phi", &GEmu_phi);
    onia_tree->Branch("GEmu_mass", &GEmu_mass);
    onia_tree->Branch("GEmu_recoId", &GEmu_recoId);

    // GEN level J/psi_1
    onia_tree->Branch("GEJpsi1_pt", &GEJpsi1_pt, "GEJpsi1_pt/D");
    onia_tree->Branch("GEJpsi1_eta", &GEJpsi1_eta, "GEJpsi1_eta/D");
    onia_tree->Branch("GEJpsi1_y", &GEJpsi1_y, "GEJpsi1_y/D");
    onia_tree->Branch("GEJpsi1_phi", &GEJpsi1_phi, "GEJpsi1_phi/D");
    onia_tree->Branch("GEJpsi1_mass", &GEJpsi1_mass, "GEJpsi1_mass/D");
    onia_tree->Branch("GEJpsi1_recoId", &GEJpsi1_recoId, "GEJpsi1_recoId/I");
    onia_tree->Branch("GEJpsi1_muId1", &GEJpsi1_muId1, "GEJpsi1_muId1/I");
    onia_tree->Branch("GEJpsi1_muId2", &GEJpsi1_muId2, "GEJpsi1_muId2/I");
    onia_tree->Branch("GEJpsi1_matchGEN", &GEJpsi1_matchGEN, "GEJpsi1_matchGEN/b");
    onia_tree->Branch("GEJpsi1_passID", &GEJpsi1_passID, "GEJpsi1_passID/b");
    onia_tree->Branch("GEJpsi1_passVtx", &GEJpsi1_passVtx, "GEJpsi1_passVtx/b");

    // GEN level J/psi_2
    onia_tree->Branch("GEJpsi2_pt", &GEJpsi2_pt, "GEJpsi2_pt/D");
    onia_tree->Branch("GEJpsi2_eta", &GEJpsi2_eta, "GEJpsi2_eta/D");
    onia_tree->Branch("GEJpsi2_y", &GEJpsi2_y, "GEJpsi2_y/D");
    onia_tree->Branch("GEJpsi2_phi", &GEJpsi2_phi, "GEJpsi2_phi/D");
    onia_tree->Branch("GEJpsi2_mass", &GEJpsi2_mass, "GEJpsi2_mass/D");
    onia_tree->Branch("GEJpsi2_recoId", &GEJpsi2_recoId, "GEJpsi2_recoId/I");
    onia_tree->Branch("GEJpsi2_muId1", &GEJpsi2_muId1, "GEJpsi2_muId1/I");
    onia_tree->Branch("GEJpsi2_muId2", &GEJpsi2_muId2, "GEJpsi2_muId2/I");
    onia_tree->Branch("GEJpsi2_matchGEN", &GEJpsi2_matchGEN, "GEJpsi2_matchGEN/b");
    onia_tree->Branch("GEJpsi2_passID", &GEJpsi2_passID, "GEJpsi2_passID/b");
    onia_tree->Branch("GEJpsi2_passVtx", &GEJpsi2_passVtx, "GEJpsi2_passVtx/b");

    // GEN level event
    onia_tree->Branch("GEevt_fourMuMass", &GEevt_fourMuMass, "GEevt_fourMuMass/D");
    onia_tree->Branch("GEevt_muPtMax", &GEevt_muPtMax, "GEevt_muPtMax/D");
    onia_tree->Branch("GEevt_valid", &GEevt_valid, "GEevt_valid/b");
    onia_tree->Branch("GEevt_passAcc", &GEevt_passAcc, "GEevt_passAcc/b");
    onia_tree->Branch("GEevt_passHLT", &GEevt_passHLT, "GEevt_passHLT/b");
    onia_tree->Branch("GEevt_matchTrg", &GEevt_matchTrg, "GEevt_matchTrg/b");
    
    // RECO level muon
    onia_tree->Branch("REmu_pt", &REmu_pt);
    onia_tree->Branch("REmu_eta", &REmu_eta);
    onia_tree->Branch("REmu_phi", &REmu_phi);
    onia_tree->Branch("REmu_mass", &REmu_mass);
    onia_tree->Branch("REmu_isSoft", &REmu_isSoft);
    onia_tree->Branch("REmu_passCut", &REmu_passCut);
    onia_tree->Branch("REmu_pvAsc", &REmu_pvAsc);

    //
    onia_tree->Branch("TEmu_pt", &TEmu_pt);
    onia_tree->Branch("TEmu_eta", &TEmu_eta);
    onia_tree->Branch("TEmu_phi", &TEmu_phi);
    onia_tree->Branch("TEmu_mass", &TEmu_mass);
    onia_tree->Branch("TEmu_pvAsc", &TEmu_pvAsc);

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
    onia_tree->Branch("REJpsi_X", &REJpsi_X);
    onia_tree->Branch("REJpsi_Y", &REJpsi_Y);
    onia_tree->Branch("REJpsi_Z", &REJpsi_Z);
    onia_tree->Branch("REJpsi_muId1", &REJpsi_muId1);
    onia_tree->Branch("REJpsi_muId2", &REJpsi_muId2);
    onia_tree->Branch("REJpsi_passCut", &REJpsi_passCut);
    // RECO level event
    onia_tree->Branch("REvtx_X", &REvtx_X);
    onia_tree->Branch("REvtx_Y", &REvtx_Y);
    onia_tree->Branch("REvtx_Z", &REvtx_Z);
    onia_tree->Branch("REevt_fourMuMass", &REevt_fourMuMass);
    onia_tree->Branch("REevt_massChisq", &REevt_massChisq);
    onia_tree->Branch("REevt_vtxProb", &REevt_vtxProb);
    onia_tree->Branch("REevt_L1muPtMax", &REevt_L1muPtMax);
    onia_tree->Branch("REevt_d", &REevt_d);
    onia_tree->Branch("REevt_X", &REevt_X);
    onia_tree->Branch("REevt_Y", &REevt_Y);
    onia_tree->Branch("REevt_Z", &REevt_Z);
    onia_tree->Branch("REevt_JpsiId1", &REevt_JpsiId1);
    onia_tree->Branch("REevt_JpsiId2", &REevt_JpsiId2);
    onia_tree->Branch("REevt_passHLT", &REevt_passHLT);
    onia_tree->Branch("REevt_matchTrg", &REevt_matchTrg);
    onia_tree->Branch("REevt_fourMuFit", &REevt_fourMuFit);
    onia_tree->Branch("REevt_twoDimuFit", &REevt_twoDimuFit);
    onia_tree->Branch("REevt_samePV", &REevt_samePV);

    prunedgenParticlesSrc_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
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
// ------------ method called for each event  ------------
void NtupleMaker::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace pat;
    Total_events_analyzed++;
	// edm::Handle<pat::CompositeCandidateCollection> dimuons;
	// iEvent.getByToken(dimuon_Label, dimuons);

	// edm::Handle<pat::CompositeCandidateCollection> conversions;
	// iEvent.getByToken(conversion_Label, conversions);

	edm::Handle<edm::View<pat::Muon>> muons;
	iEvent.getByToken(muon_Label, muons);

    edm::Handle<edm::View<pat::PackedCandidate>> candidates;
    iEvent.getByToken(candidate_Label, candidates);

	edm::Handle<reco::VertexCollection> primaryVertices_handle;
	iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);
	reco::Vertex thePrimaryV;

	edm::Handle<reco::BeamSpot> theBeamSpot;
	iEvent.getByToken(bs_Label,theBeamSpot);
	reco::BeamSpot bs = *theBeamSpot;

	reco::VertexCollection primaryVertices = *primaryVertices_handle;
    primaryVertices.push_back(reco::Vertex(bs.position(), bs.covariance3D()));
    REvtx_X.clear();
    REvtx_Y.clear();
    REvtx_Z.clear();
    for(auto it = primaryVertices.begin(); it != primaryVertices.end(); ++it) {
        REvtx_X.push_back(it->x());
        REvtx_Y.push_back(it->y());
        REvtx_Z.push_back(it->z());
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

    edm::Handle<edm::TriggerResults> hltR;
    iEvent.getByToken(triggerResultsTok_, hltR);
    edm::TriggerNames trgNames = iEvent.triggerNames(*hltR);
    int nTrig = hltR->size();
    evtPassHLT = false;
    for(int i = 0; i < nTrig; i++) {// for each of existing trigger path
        string trgName = trgNames.triggerName(i);
        if(!(*hltR)[i].wasrun() || !(*hltR)[i].accept() || (*hltR)[i].error()) continue;
        bool foundTrg = false;
        for(int j = 0; j < (int)triggerList.size(); j++) {// if any of wanted trigger paths included
            if(strstr(trgName.c_str(), triggerList[j].c_str())) {
                foundTrg = true;
                evtPassHLT = true;
                Total_events_triggered++;
                break;
            }
        }
        // vout<<"Investigating: "<<trgName<<endl;
        if(foundTrg) break;
    }

    run     = iEvent.id().run();
    lumi    = iEvent.id().luminosityBlock();
    event   = iEvent.id().event();
    vout<<"run: "<<run<<endl;
    vout<<"lumi: "<<lumi<<endl;
    vout<<"event: "<<event<<endl;

    GEmu_pt.clear();
    GEmu_eta.clear();
    GEmu_phi.clear();
    GEmu_mass.clear();
    GEmu_recoId.clear();

    GEJpsi1_recoId = -1;
    GEJpsi1_passVtx = false;
    GEJpsi2_recoId = -1;
    GEJpsi2_passVtx = false;

    GEevt_fourMuMass = 0;
    GEevt_muPtMax = 0;
    GEevt_passHLT = false;
    GEevt_matchTrg = false;
    
    REmu_pt.clear();
    REmu_eta.clear();
    REmu_phi.clear();
    REmu_mass.clear();
    REmu_isSoft.clear();
    REmu_passCut.clear();
    REmu_pvAsc.clear();

    TEmu_pt.clear();
    TEmu_eta.clear();
    TEmu_phi.clear();
    TEmu_mass.clear();
    TEmu_pvAsc.clear();

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
    REJpsi_X.clear();
    REJpsi_Y.clear();
    REJpsi_Z.clear();
    REJpsi_muId1.clear();
    REJpsi_muId2.clear();
    REJpsi_passCut.clear();
    
    REevt_fourMuMass.clear();
    REevt_massChisq.clear();
    REevt_vtxProb.clear();
    REevt_L1muPtMax.clear();
    REevt_d.clear();
    REevt_X.clear();
    REevt_Y.clear();
    REevt_Z.clear();
    REevt_JpsiId1.clear();
    REevt_JpsiId2.clear();
    REevt_passHLT.clear();
    REevt_matchTrg.clear();
    REevt_fourMuFit.clear();
    REevt_twoDimuFit.clear();
    REevt_samePV.clear();

	edm::ESHandle<MagneticField> bFieldHandle;
	iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

    edm::Handle<reco::GenParticleCollection> prunedgenParticles;
    iEvent.getByToken(prunedgenParticlesSrc_, prunedgenParticles);

	if(prunedgenParticles.isValid()) {
        vout<<"setting gen variables"<<endl;    
        setGENVariables(prunedgenParticles);
	}

    srand(time(0));
    fourMuonFit(iEvent, hltR, muons, candidates, bFieldHandle, bs, thePrimaryV,primaryVertices_handle);
    
    onia_tree->Fill();
}
// ------------ method called once each job just before starting event loop  ------------
void NtupleMaker::beginJob() {
    vout<<"The job Begins"<<endl;
}
// ------------ method called once each job just after ending the event loop  ------------
void NtupleMaker::endJob() {
	cout <<"Job End"<<endl;
	cout<<"Total events analyzed: "<<Total_events_analyzed<<endl;
	cout<<"Total events triggered: "<<Total_events_triggered<<endl;
	cout<<"Total events trg-matched: "<<Total_events_trg_matched<<endl;
}
// ------------ method called when ending the processing of a run  ------------
void NtupleMaker::beginRun(edm::Run & iRun, edm::EventSetup const& iSetup) {
    vout<<"New run"<<endl;
    //--- m_l1GtUtils.getL1GtRunCache(run, iSetup, true, false);
    bool hltConfigChanged, test = hltConfig_.init(iRun, iSetup, hltName_, hltConfigChanged);
    if(test) {
        vout<<"hltConfig_.size() = "<<hltConfig_.size()<<endl;
        // check iftrigger name in (new) config
        const unsigned int n(hltConfig_.size());
        const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
        vout<<" triggerIndex = "<<triggerIndex<<endl;
        if(triggerIndex >= n) {
            vout<<"HLT event analyzer MINIAOD beginRun:"<<" TriggerName "<<triggerName_<<" not available in (new) config!"<<endl;
            vout<<"Available TriggerNames are: "<<endl;
            hltConfig_.dump("Triggers");
        }
    }else vout<<"beginRun:"<<" HLT config extraction failure with process name "<<hltName_<<endl;
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

bool NtupleMaker::triggerMatch(const edm::Event &iEvent, const edm::Handle<edm::TriggerResults>& hltR, pat::Muon& muon, enum TrgLabel flag) {
    const pat::TriggerObjectStandAloneCollection& trgObjList = muon.triggerObjectMatches();
    vout<<"Trigger matching:"<<(int)trgObjList.size()<<endl;
    for(int i = 0; i < (int)trgObjList.size(); i++) {
        pat::TriggerObjectStandAlone trgObj(trgObjList[i]);
        trgObj.unpackNamesAndLabels(iEvent, *hltR);
        const vector<string>& fltLabels = trgObj.filterLabels();
        for(int j = 0; j < (int)fltLabels.size() && j < 10; j++) vout<<"Label: "<<fltLabels[j]<<endl;
        if(flag == Vtx) {
            if(trgObj.filter("hltVertexmumuFilterJpsiMuon")) return true;
        }else if(flag == L3) {
            if(trgObj.filter("hltTripleMuL3PreFiltered0")) return true;
        }else if(flag == L1) {
            if(trgObj.filter("hltL1sTripleMu0")) return true;
        }
    }
    return false;
}
void NtupleMaker::setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles) {
    vector<const reco::Candidate*> temp_Jpsi, temp_JpsiMom;
    vector<vector<const reco::Candidate*>> temp_JpsiDau;
    bool passAcc = true, valid = false;
    TLorentzVector temp_fourMuLV;
    for(auto p = prunedgenParticles->begin(); p != prunedgenParticles->end(); p++) {
        if(p->pdgId() != 443 || p->status() != 2) continue;
        vector<const reco::Candidate*> temp_mu;
        for(int i = 0; i < (int)p->numberOfDaughters(); i++) {
            const reco::Candidate *p_dau = p->daughter(i);
            while(p_dau->status() != 1 && p_dau->daughter(0)) p_dau = p_dau->daughter(0);
            if(fabs(p_dau->pdgId()) == 13 && p_dau->status() == 1) temp_mu.push_back(p_dau);//if daughter of J/psi is muon
        }
        if((int)temp_mu.size() != 2) continue;
        temp_JpsiDau.push_back(temp_mu);
        const reco::Candidate *Jpsi_self = &*p;
        temp_Jpsi.push_back(Jpsi_self);
        while(Jpsi_self->mother(0)->pdgId() == 443) Jpsi_self = Jpsi_self->mother(0);
        temp_JpsiMom.push_back(Jpsi_self->mother(0));
    }
    for(int i = 0; i < (int)temp_Jpsi.size(); i++) {
        for(int j = i + 1; j < (int)temp_Jpsi.size(); j++) {
            if(temp_JpsiMom[i] != temp_JpsiMom[j]) continue;//only for SPS
            valid = true;
            for(int k = 0; k < 2; k++) {
                GEmu_pt.push_back(temp_JpsiDau[i][k]->pt());
                GEmu_eta.push_back(temp_JpsiDau[i][k]->eta());
                GEmu_phi.push_back(temp_JpsiDau[i][k]->phi());
                GEmu_mass.push_back(temp_JpsiDau[i][k]->mass());
            }
            for(int k = 0; k < 2; k++) {
                GEmu_pt.push_back(temp_JpsiDau[j][k]->pt());
                GEmu_eta.push_back(temp_JpsiDau[j][k]->eta());
                GEmu_phi.push_back(temp_JpsiDau[j][k]->phi());
                GEmu_mass.push_back(temp_JpsiDau[j][k]->mass());
            }
            GEJpsi1_pt = temp_Jpsi[i]->pt();
            GEJpsi1_eta = temp_Jpsi[i]->eta();
            GEJpsi1_y = temp_Jpsi[i]->rapidity();
            GEJpsi1_phi = temp_Jpsi[i]->phi();
            GEJpsi1_mass = temp_Jpsi[i]->mass();
            GEJpsi1_muId1 = 0;
            GEJpsi1_muId2 = 1;
            GEJpsi2_pt = temp_Jpsi[j]->pt();
            GEJpsi2_eta = temp_Jpsi[j]->eta();
            GEJpsi2_y = temp_Jpsi[j]->rapidity();
            GEJpsi2_phi = temp_Jpsi[j]->phi();
            GEJpsi2_mass = temp_Jpsi[j]->mass();
            GEJpsi2_muId1 = 2;
            GEJpsi2_muId2 = 3;
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
    const edm::Event &iEvent,
    const edm::Handle<edm::TriggerResults>& hltR,
    edm::Handle<edm::View<pat::Muon>> muons,
    edm::Handle<edm::View<pat::PackedCandidate>> candidates,
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
    if(nSelMuon < 4) return;
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
    GEJpsi1_matchGEN = GEevt_valid && GEmu_recoId[GEJpsi1_muId1] != -1 && GEmu_recoId[GEJpsi1_muId2] != -1;
    GEJpsi2_matchGEN = GEevt_valid && GEmu_recoId[GEJpsi2_muId1] != -1 && GEmu_recoId[GEJpsi2_muId2] != -1;
    // single muon cut
    GEJpsi1_passID = GEJpsi1_matchGEN && REmu_passCut[GEmu_recoId[GEJpsi1_muId1]] && REmu_passCut[GEmu_recoId[GEJpsi1_muId2]];
    GEJpsi2_passID = GEJpsi2_matchGEN && REmu_passCut[GEmu_recoId[GEJpsi2_muId1]] && REmu_passCut[GEmu_recoId[GEJpsi2_muId2]];
    // relate muons to packed candidates and primary vertices
    vector<const pat::PackedCandidate*> canMuons;
    for(edm::View<pat::PackedCandidate>::const_iterator recoCand = candidates->begin(); recoCand != candidates->end(); ++recoCand) {
        if(!recoCand->isGlobalMuon()) continue;
        canMuons.push_back(&(*recoCand));
        TEmu_pt.push_back(recoCand->pt());
        TEmu_eta.push_back(recoCand->eta());
        TEmu_phi.push_back(recoCand->phi());
        TEmu_mass.push_back(recoCand->mass());
        TEmu_pvAsc.push_back((int)(recoCand->vertexRef().key()));
    }
    for(int i = 0; i < nSelMuon; i++) {
        int id = -1;
        double min = 0.1;
        for(int j = 0; j < (int)canMuons.size(); j++) {
            if(fabs(REmu_pt[i] - canMuons[j]->pt()) > min) continue;
            min = fabs(REmu_pt[i] - canMuons[j]->pt());
            id = j;
        }
        if(id == -1) {
            REmu_pvAsc.push_back(-1);// represents no associated candidate
            continue;
        }
        REmu_pvAsc.push_back((int)(canMuons[id]->vertexRef().key()));
        canMuons.erase(canMuons.begin() + id);
    }
    // combinition of muons
    KinematicConstraint *JpsiMassCons = new MassKinematicConstraint(JpsiMass, JpsiSigma);
    vector<RefCountedKinematicTree> JpsiVtxTrees;// Save J/psi,psi(2S) vertex tree
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
            if(selDimuonLV.M() < 2.7 || selDimuonLV.M() > 3.5) continue;
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
            bool dimuCutPartial = dimuonVtxProb > 0.005 && fabs(temp_dimuLV.Pt() - 25) < 15 && fabs(temp_dimuLV.Rapidity()) < 2;//temp_dimuLV.Pt() > 5
            // dimuon vertex fit with mass constraint
            if(dimuonMass < 2.7 || dimuonMass > 3.5) continue;// Current dimuon matches J/psi
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
            REJpsi_X.push_back(vFitVertex->position().x());
            REJpsi_Y.push_back(vFitVertex->position().y());
            REJpsi_Z.push_back(vFitVertex->position().z());
            REJpsi_muId1.push_back(muId1);
            REJpsi_muId2.push_back(muId2);
            REJpsi_passCut.push_back(dimuCutPartial && dimuonMass > 2.95 && dimuonMass < 3.25);
            if(GEJpsi1_passID && ((GEmu_recoId[GEJpsi1_muId1] == i && GEmu_recoId[GEJpsi1_muId2] == j) || (GEmu_recoId[GEJpsi1_muId1] == j && GEmu_recoId[GEJpsi1_muId2] == i))) {
                GEJpsi1_recoId = REJpsi_pt.size() - 1;
                GEJpsi1_passVtx = REJpsi_passCut[GEJpsi1_recoId];
            }
            if(GEJpsi2_passID && ((GEmu_recoId[GEJpsi2_muId1] == i && GEmu_recoId[GEJpsi2_muId2] == j) || (GEmu_recoId[GEJpsi2_muId1] == j && GEmu_recoId[GEJpsi2_muId2] == i))) {
                GEJpsi2_recoId = REJpsi_pt.size() - 1;
                GEJpsi2_passVtx = REJpsi_passCut[GEJpsi2_recoId];
            }
            JpsiVtxTrees.push_back(dimuonVtxFitTree);// Save J/psi vertex tree
        }
    }
    int nJpsi = REJpsi_pt.size();
    vout<<"In this event: ["<<nJpsi<<"] J/psi"<<endl;
    if(nJpsi < 2) return;
    // check HLT
    if (evtPassHLT) vout<<"Event passed trigger selection"<<endl;
    // else return;
    GEevt_passHLT = GEJpsi1_passVtx && GEJpsi2_passVtx && evtPassHLT;
    // four muon vertex fit and trigger matching
    for(int i = 0; i < nJpsi; i++) {
        if(!REJpsi_passCut[i]) continue;
        int i1 = REJpsi_muId1[i], i2 = REJpsi_muId2[i];
        for(int j = i + 1; j < nJpsi; j++) {
            if(!REJpsi_passCut[j]) continue;
            int j1 = REJpsi_muId1[j], j2 = REJpsi_muId2[j];
            if(i1 == j1 || i1 == j2 || i2 == j1 || i2 == j2) continue;
            RefCountedKinematicTree JpsiVtxTree1 = JpsiVtxTrees[i], JpsiVtxTree2 = JpsiVtxTrees[j];
            JpsiVtxTree1->movePointerToTheTop();
            JpsiVtxTree2->movePointerToTheTop();
            // trigger matching
            int nVtxMuon = triggerMatch(iEvent, hltR, selMuons[i1], Vtx) + triggerMatch(iEvent, hltR, selMuons[i2], Vtx);
            int nL3Muon = triggerMatch(iEvent, hltR, selMuons[j1], L3) + triggerMatch(iEvent, hltR, selMuons[j2], L3);
            bool evtTrgMatch = nVtxMuon >= 2 && nL3Muon >= 1;
            nVtxMuon = triggerMatch(iEvent, hltR, selMuons[j1], Vtx) + triggerMatch(iEvent, hltR, selMuons[j2], Vtx);
            nL3Muon = triggerMatch(iEvent, hltR, selMuons[i1], L3) + triggerMatch(iEvent, hltR, selMuons[i2], L3);
            evtTrgMatch = evtTrgMatch || (nVtxMuon >= 2 && nL3Muon >= 1);
            if(evtTrgMatch) {
                Total_events_trg_matched++;
            }
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
            Double_t fourMuMass = 0, fourMuVtxProb = -1, fourMuX = 0, fourMuY = 0, fourMuZ = 0;
            bool passFitFourMu = !fourMuVtxFitTree->isEmpty();
            if(passFitFourMu) {
                fourMuVtxFitTree->movePointerToTheTop();
                RefCountedKinematicParticle fourMu = fourMuVtxFitTree->currentParticle();
                if(fourMu->currentState().isValid()) fourMuMass = fourMu->currentState().mass();
                RefCountedKinematicVertex fourMuVtx = fourMuVtxFitTree->currentDecayVertex();
                if(fourMuVtx->vertexIsValid()) {
                    fourMuVtxProb = ChiSquaredProbability((double)(fourMuVtx->chiSquared()), (double)(fourMuVtx->degreesOfFreedom()));
                    fourMuX = fourMuVtx->position().x();
                    fourMuY = fourMuVtx->position().y();
                    fourMuZ = fourMuVtx->position().z();
                }
            }
            vout<<"Four-muon vertex probability: "<<fourMuVtxProb<<endl;
            // Method 2: split into 2 dimuons
            vector<RefCountedKinematicParticle> fitDimus;
            RefCountedKinematicParticle JpsiVtx1 = JpsiVtxTree1->currentParticle(), JpsiVtx2 = JpsiVtxTree2->currentParticle();
            fitDimus.push_back(JpsiVtx1);
            fitDimus.push_back(JpsiVtx2);
            KinematicParticleVertexFitter twoDimuFitter;
            RefCountedKinematicTree twoDimuVtxFitTree = twoDimuFitter.fit(fitDimus);
            bool passFitDimu = !twoDimuVtxFitTree->isEmpty();
            // Handel output variables
            vector<Double_t> temp_L1muPt;
            if(triggerMatch(iEvent, hltR, selMuons[i1], L1)) temp_L1muPt.push_back(REmu_pt[i1]);
            if(triggerMatch(iEvent, hltR, selMuons[i2], L1)) temp_L1muPt.push_back(REmu_pt[i2]);
            if(triggerMatch(iEvent, hltR, selMuons[j1], L1)) temp_L1muPt.push_back(REmu_pt[j1]);
            if(triggerMatch(iEvent, hltR, selMuons[j2], L1)) temp_L1muPt.push_back(REmu_pt[j2]);
            Double_t temp_L1muPtMax = temp_L1muPt.empty() ? 0 : *max_element(temp_L1muPt.begin(), temp_L1muPt.end());
            Double_t d = calD(JpsiVtxTree1->currentDecayVertex(), JpsiVtxTree2->currentDecayVertex());
            Double_t massChisq = SQUARE((REJpsi_mass[i] - JpsiMass) / REJpsi_massErr[i]) + SQUARE((REJpsi_mass[j] - JpsiMass) / REJpsi_massErr[j]);
            bool samePV = (REmu_pvAsc[i1] == REmu_pvAsc[i2] && REmu_pvAsc[i2] == REmu_pvAsc[j1] && REmu_pvAsc[j1] == REmu_pvAsc[j2]);
            int pos = lower_bound(REevt_massChisq.begin(), REevt_massChisq.end(), massChisq, greater<Double_t>()) - REevt_massChisq.begin();
            REevt_fourMuMass.insert(REevt_fourMuMass.begin() + pos, fourMuMass);
            REevt_massChisq.insert(REevt_massChisq.begin() + pos, massChisq);
            REevt_vtxProb.insert(REevt_vtxProb.begin() + pos, fourMuVtxProb);
            REevt_L1muPtMax.insert(REevt_L1muPtMax.begin() + pos, temp_L1muPtMax);
            REevt_d.insert(REevt_d.begin() + pos, d);
            REevt_X.insert(REevt_X.begin() + pos, fourMuX);
            REevt_Y.insert(REevt_Y.begin() + pos, fourMuY);
            REevt_Z.insert(REevt_Z.begin() + pos, fourMuZ);
            if((double)rand() / (RAND_MAX + 1.0) > 0.5) {
                REevt_JpsiId1.insert(REevt_JpsiId1.begin() + pos, i);
                REevt_JpsiId2.insert(REevt_JpsiId2.begin() + pos, j);
            }else {
                REevt_JpsiId1.insert(REevt_JpsiId1.begin() + pos, j);
                REevt_JpsiId2.insert(REevt_JpsiId2.begin() + pos, i);
            }
            REevt_passHLT.insert(REevt_passHLT.begin() + pos, evtPassHLT);
            REevt_matchTrg.insert(REevt_matchTrg.begin() + pos, evtTrgMatch);
            REevt_fourMuFit.insert(REevt_fourMuFit.begin() + pos, passFitFourMu);
            REevt_twoDimuFit.insert(REevt_twoDimuFit.begin() + pos, passFitDimu);
            REevt_samePV.insert(REevt_samePV.begin() + pos, samePV);
            if((GEJpsi1_recoId == i && GEJpsi2_recoId == j) || (GEJpsi1_recoId == j && GEJpsi2_recoId == i)) {
                GEevt_matchTrg = GEevt_passHLT && evtTrgMatch && samePV;
                GEevt_muPtMax = max(REmu_pt[i1], REmu_pt[j1]);
            }
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
