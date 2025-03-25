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

#define VERBOSE// switch between verbose output and concise output
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
const float psi2SMass = 3.686097;
const float JpsiSigma = 0.00004;
const float psi2SSigma = 0.00004;
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
        void fourMuonFit(const edm::Event &iEvent, const edm::Handle<edm::TriggerResults>& hltR, edm::Handle< edm::View<pat::Muon> > muons, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV, edm::Handle<reco::VertexCollection>);
        void calCtau(RefCountedKinematicVertex& decayVrtx, RefCountedKinematicParticle& kinePart, Vertex& bs, Double_t& ctau, Double_t& ctauErr);

        virtual void beginJob();
        virtual void analyze(const edm::Event &, const edm::EventSetup &);
        virtual void endJob();

        virtual void beginRun(edm::Run &, edm::EventSetup const&);
        virtual void endRun(edm::Run const &, edm::EventSetup const &);
        virtual void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);
        virtual void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);
		// ----------member data ---------------------------
		edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
		edm::EDGetTokenT<reco::BeamSpot> bs_Label;
		edm::EDGetTokenT<edm::View<pat::Muon>> muon_Label;
		edm::EDGetTokenT<edm::TriggerResults> triggerResultsTok_;
        edm::EDGetTokenT<trigger::TriggerEvent> triggerEventTok_;

        // RECO level muon
        vector<Double_t> REmu_pt, REmu_eta, REmu_phi, REmu_mass;
        vector<bool> REmu_isSoft, REmu_passCut;
        // RECO level J/psi
        vector<Double_t> REJpsi_pt, REJpsi_eta, REJpsi_phi, REJpsi_mass, REJpsi_massErr, REJpsi_vtxProb, REJpsi_cstrVtxProb, REJpsi_ctau, REJpsi_ctauErr;
        vector<int> REJpsi_muId1, REJpsi_muId2;
        vector<bool> REJpsi_passCut;
        // RECO level psi(2S)
        vector<Double_t> REpsi2S_pt, REpsi2S_eta, REpsi2S_phi, REpsi2S_mass, REpsi2S_massErr, REpsi2S_vtxProb, REpsi2S_cstrVtxProb, REpsi2S_ctau, REpsi2S_ctauErr;
        vector<int> REpsi2S_muId1, REpsi2S_muId2;
        vector<bool> REpsi2S_passCut;
        // RECO level event
        vector<Double_t> REevt_fourMuMass, REevt_massChisq, REevt_vtxProb, REevt_L1muPtMax;
        vector<int> REevt_JpsiId, REevt_psi2SId;
        vector<bool> REevt_passHLT, REevt_matchTrg;

        bool evtPassHLT, hasMatchEvt;
        // some counters for cut flow
        int Total_events_analyzed;
        int Total_events_triggered;
        int Total_events_trg_matched;
        
		ULong64_t run;
		ULong64_t lumi;
        ULong64_t event;
        vector<Float_t> primaryVtx;
 
        string hltName_;
        string triggerName_;
        vector<string> triggerList;
        HLTConfigProvider hltConfig_;
		TTree *onia_tree;
};
//
// constructors and destructor
//
NtupleMaker::NtupleMaker(const edm::ParameterSet & iConfig):
    primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
    bs_Label(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("offlineBeamSpot"))),
    muon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter< edm::InputTag>("muons"))),
    triggerResultsTok_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
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
    onia_tree->Branch("REJpsi_phi", &REJpsi_phi);
    onia_tree->Branch("REJpsi_mass", &REJpsi_mass);
    onia_tree->Branch("REJpsi_massErr", &REJpsi_massErr);
    onia_tree->Branch("REJpsi_vtxProb", &REJpsi_vtxProb);
    onia_tree->Branch("REJpsi_cstrVtxProb", &REJpsi_cstrVtxProb);
    onia_tree->Branch("REJpsi_ctau", &REJpsi_ctau);
    onia_tree->Branch("REJpsi_ctauErr", &REJpsi_ctauErr);
    onia_tree->Branch("REJpsi_muId1", &REJpsi_muId1);
    onia_tree->Branch("REJpsi_muId2", &REJpsi_muId2);
    onia_tree->Branch("REJpsi_passCut", &REJpsi_passCut);
    // RECO level psi(2S)
    onia_tree->Branch("REpsi2S_pt", &REpsi2S_pt);
    onia_tree->Branch("REpsi2S_eta", &REpsi2S_eta);
    onia_tree->Branch("REpsi2S_phi", &REpsi2S_phi);
    onia_tree->Branch("REpsi2S_mass", &REpsi2S_mass);
    onia_tree->Branch("REpsi2S_massErr", &REpsi2S_massErr);
    onia_tree->Branch("REpsi2S_vtxProb", &REpsi2S_vtxProb);
    onia_tree->Branch("REpsi2S_cstrVtxProb", &REpsi2S_cstrVtxProb);
    onia_tree->Branch("REpsi2S_ctau", &REpsi2S_ctau);
    onia_tree->Branch("REpsi2S_ctauErr", &REpsi2S_ctauErr);
    onia_tree->Branch("REpsi2S_muId1", &REpsi2S_muId1);
    onia_tree->Branch("REpsi2S_muId2", &REpsi2S_muId2);
    onia_tree->Branch("REpsi2S_passCut", &REpsi2S_passCut);
    // RECO level event
    onia_tree->Branch("REevt_fourMuMass", &REevt_fourMuMass);
    onia_tree->Branch("REevt_massChisq", &REevt_massChisq);
    onia_tree->Branch("REevt_vtxProb", &REevt_vtxProb);
    onia_tree->Branch("REevt_L1muPtMax", &REevt_L1muPtMax);
    onia_tree->Branch("REevt_JpsiId", &REevt_JpsiId);
    onia_tree->Branch("REevt_psi2SId", &REevt_psi2SId);
    onia_tree->Branch("REevt_passHLT", &REevt_passHLT);
    onia_tree->Branch("REevt_matchTrg", &REevt_matchTrg);

    onia_tree->Branch("primaryVtx", &primaryVtx);
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

	edm::Handle<edm::View<pat::Muon>> muons;
	iEvent.getByToken(muon_Label, muons);

	edm::Handle<reco::VertexCollection> primaryVertices_handle;
	iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);
	reco::Vertex thePrimaryV;

	edm::Handle<reco::BeamSpot> theBeamSpot;
	iEvent.getByToken(bs_Label,theBeamSpot);
	reco::BeamSpot bs = *theBeamSpot;

	if(primaryVertices_handle->begin() != primaryVertices_handle->end()) thePrimaryV = reco::Vertex(*(primaryVertices_handle->begin()));
	else thePrimaryV = reco::Vertex(bs.position(), bs.covariance3D());
	primaryVtx.clear();
    primaryVtx.push_back(thePrimaryV.x());
    primaryVtx.push_back(thePrimaryV.y());
    primaryVtx.push_back(thePrimaryV.z());

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
    
    REmu_pt.clear();
    REmu_eta.clear();
    REmu_phi.clear();
    REmu_mass.clear();
    REmu_isSoft.clear();
    REmu_passCut.clear();
    REJpsi_pt.clear();
    REJpsi_eta.clear();
    REJpsi_phi.clear();
    REJpsi_mass.clear();
    REJpsi_massErr.clear();
    REJpsi_vtxProb.clear();
    REJpsi_ctau.clear();
    REJpsi_ctauErr.clear();
    REJpsi_muId1.clear();
    REJpsi_muId2.clear();
    REJpsi_passCut.clear();
    REpsi2S_pt.clear();
    REpsi2S_eta.clear();
    REpsi2S_phi.clear();
    REpsi2S_mass.clear();
    REpsi2S_massErr.clear();
    REpsi2S_vtxProb.clear();
    REpsi2S_cstrVtxProb.clear();
    REpsi2S_ctau.clear();
    REpsi2S_ctauErr.clear();
    REpsi2S_muId1.clear();
    REpsi2S_muId2.clear();
    REpsi2S_passCut.clear();
    REevt_fourMuMass.clear();
    REevt_massChisq.clear();
    REevt_vtxProb.clear();
    REevt_L1muPtMax.clear();
    REevt_JpsiId.clear();
    REevt_psi2SId.clear();
    REevt_passHLT.clear();
    REevt_matchTrg.clear();

    hasMatchEvt = false;

	edm::ESHandle<MagneticField> bFieldHandle;
	iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
    RefCountedKinematicTree tree_ups_part1;

    fourMuonFit(iEvent, hltR, muons, bFieldHandle, bs, thePrimaryV,primaryVertices_handle);
    
    if(hasMatchEvt) onia_tree->Fill();
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
            if(trgObj.filter("hltVertexmumuFilterJpsiMuon3p5")) return true;
        }else if(flag == L3) {
            if(trgObj.filter("hltTripleMuL3PreFiltered222")) return true;
        }else if(flag == L1) {
            if(trgObj.filter("hltL1TripleMu5SQ3SQ0OQDoubleMu53SQOSMassMax9")) return true;
        }
    }
    return false;
}

void NtupleMaker::fourMuonFit(
    const edm::Event &iEvent,
    const edm::Handle<edm::TriggerResults>& hltR,
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
    if(nSelMuon < 4) return;
    // combinition of muons
    KinematicConstraint *JpsiMassCons = new MassKinematicConstraint(JpsiMass, JpsiSigma), *psi2SMassCons = new MassKinematicConstraint(psi2SMass, psi2SSigma);
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
            bool matchJpsi = dimuonMass > 2.7 && dimuonMass < 3.5;
            bool matchpsi2S = dimuonMass > 3.3 && dimuonMass < 4.1;
            TLorentzVector temp_dimuLV;
            temp_dimuLV.SetPxPyPzE(vFitDimuState.globalMomentum().x(), vFitDimuState.globalMomentum().y(), vFitDimuState.globalMomentum().z(), vFitDimuState.kinematicParameters().energy());
            Double_t ctau, ctauErr;
            calCtau(vFitVertex, vFitDimu, thePrimaryV, ctau, ctauErr);
            int muId1 = i, muId2 = i;
            if(REmu_pt[i] > REmu_pt[j]) muId2 = j;
            else muId1 = j;
            bool dimuCutPartial = dimuonVtxProb > 0.005 && temp_dimuLV.Pt() > 5 && fabs(temp_dimuLV.Eta()) < 2;
            // dimuon vertex fit with mass constraint
            if(matchJpsi) {
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
                REJpsi_phi.push_back(temp_dimuLV.Phi());
                REJpsi_mass.push_back(dimuonMass);
                REJpsi_massErr.push_back(dimuonMassErr);
                REJpsi_vtxProb.push_back(dimuonVtxProb);
                REJpsi_cstrVtxProb.push_back(JpsiCstrVtxProb);
                REJpsi_ctau.push_back(ctau);
                REJpsi_ctauErr.push_back(ctauErr);
                REJpsi_muId1.push_back(muId1);
                REJpsi_muId2.push_back(muId2);
                REJpsi_passCut.push_back(dimuCutPartial && dimuonMass > 2.85 && dimuonMass < 3.35);
            }
            if(matchpsi2S) {
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
                REpsi2S_phi.push_back(temp_dimuLV.Phi());
                REpsi2S_mass.push_back(dimuonMass);
                REpsi2S_massErr.push_back(dimuonMassErr);
                REpsi2S_vtxProb.push_back(dimuonVtxProb);
                REpsi2S_cstrVtxProb.push_back(psi2SCstrVtxProb);
                REpsi2S_ctau.push_back(ctau);
                REpsi2S_ctauErr.push_back(ctauErr);
                REpsi2S_muId1.push_back(muId1);
                REpsi2S_muId2.push_back(muId2);
                REpsi2S_passCut.push_back(dimuCutPartial && dimuonMass > 3.35 && dimuonMass < 4.05);
            }
        }
    }
    int nJpsi = REJpsi_pt.size(), npsi2S = REpsi2S_pt.size();
    vout<<"In this event: ["<<nJpsi<<"] J/psi, ["<<npsi2S<<"] psi(2S)"<<endl;
    if(!nJpsi || !npsi2S) return;
    // check HLT
    if (evtPassHLT) {vout<<"Event passed trigger selection"<<endl;}
    else return;
    // four muon vertex fit and trigger matching
    for(int i = 0; i < nJpsi; i++) {
        if(!REJpsi_passCut[i]) continue;
        int i1 = REJpsi_muId1[i], i2 = REJpsi_muId2[i];
        for(int j = 0; j < npsi2S; j++) {
            if(!REpsi2S_passCut[j]) continue;
            int j1 = REpsi2S_muId1[j], j2 = REpsi2S_muId2[j];
            if(i1 == j1 || i1 == j2 || i2 == j1 || i2 == j2) continue;
            // four muon vertex fit
            reco::TransientTrack muoni1TrsTrack(selMuons[i1].track(), &(*bFieldHandle));
            reco::TransientTrack muoni2TrsTrack(selMuons[i2].track(), &(*bFieldHandle));
            reco::TransientTrack muonj1TrsTrack(selMuons[j1].track(), &(*bFieldHandle));
            reco::TransientTrack muonj2TrsTrack(selMuons[j2].track(), &(*bFieldHandle));
            KinematicParticleFactoryFromTransientTrack factory;
            vector<RefCountedKinematicParticle> fitMuons;
            fitMuons.push_back(factory.particle(muoni1TrsTrack, muonMass, 0., 0., muonSigma));
            fitMuons.push_back(factory.particle(muoni2TrsTrack, muonMass, 0., 0., muonSigma));
            fitMuons.push_back(factory.particle(muonj1TrsTrack, muonMass, 0., 0., muonSigma));
            fitMuons.push_back(factory.particle(muonj2TrsTrack, muonMass, 0., 0., muonSigma));
            KinematicParticleVertexFitter fourMuFitter;
            RefCountedKinematicTree fourMuVtxFitTree;
            fourMuVtxFitTree = fourMuFitter.fit(fitMuons);
            if(fourMuVtxFitTree->isEmpty()) continue;
            fourMuVtxFitTree->movePointerToTheTop();
            RefCountedKinematicParticle fourMu = fourMuVtxFitTree->currentParticle();
            RefCountedKinematicVertex fourMuVtx = fourMuVtxFitTree->currentDecayVertex();
            Double_t fourMuMass = fourMu->currentState().mass();
            Double_t fourMuVtxProb = ChiSquaredProbability((double)(fourMuVtx->chiSquared()), (double)(fourMuVtx->degreesOfFreedom()));
            vout<<"Four-muon vertex probability: "<<fourMuVtxProb<<endl;
            if(fourMuMass < 7.5) continue;
            // trigger matching
            int nVtxMuon = triggerMatch(iEvent, hltR, selMuons[i1], Vtx) + triggerMatch(iEvent, hltR, selMuons[i2], Vtx);
            int nL3Muon = triggerMatch(iEvent, hltR, selMuons[j1], L3) + triggerMatch(iEvent, hltR, selMuons[j2], L3);
            bool evtTrgMatch = nVtxMuon >= 2 && nL3Muon >= 1;
            if(evtTrgMatch) {
                Total_events_trg_matched++;
                hasMatchEvt = true;
            }
            vector<Double_t> temp_L1muPt;
            if(triggerMatch(iEvent, hltR, selMuons[i1], L1)) temp_L1muPt.push_back(REmu_pt[i1]);
            if(triggerMatch(iEvent, hltR, selMuons[i2], L1)) temp_L1muPt.push_back(REmu_pt[i2]);
            if(triggerMatch(iEvent, hltR, selMuons[j1], L1)) temp_L1muPt.push_back(REmu_pt[j1]);
            if(triggerMatch(iEvent, hltR, selMuons[j2], L1)) temp_L1muPt.push_back(REmu_pt[j2]);
            Double_t temp_L1muPtMax = temp_L1muPt.empty() ? 0 : *max_element(temp_L1muPt.begin(), temp_L1muPt.end());
            Double_t massChisq = SQUARE((REJpsi_mass[i] - JpsiMass) / REJpsi_massErr[i]) + SQUARE((REpsi2S_mass[j] - psi2SMass) / REpsi2S_massErr[j]);
            int pos = lower_bound(REevt_massChisq.begin(), REevt_massChisq.end(), massChisq, greater<Double_t>()) - REevt_massChisq.begin();
            REevt_fourMuMass.insert(REevt_fourMuMass.begin() + pos, fourMuMass);
            REevt_massChisq.insert(REevt_massChisq.begin() + pos, massChisq);
            REevt_vtxProb.insert(REevt_vtxProb.begin() + pos, fourMuVtxProb);
            REevt_L1muPtMax.insert(REevt_L1muPtMax.begin() + pos, temp_L1muPtMax);
            REevt_JpsiId.insert(REevt_JpsiId.begin() + pos, i);
            REevt_psi2SId.insert(REevt_psi2SId.begin() + pos, j);
            REevt_passHLT.insert(REevt_passHLT.begin() + pos, evtPassHLT);
            REevt_matchTrg.insert(REevt_matchTrg.begin() + pos, evtTrgMatch);
        }
    }
    return;
}

void NtupleMaker::calCtau(RefCountedKinematicVertex& decayVrtx, RefCountedKinematicParticle& kinePart, Vertex& bs, Double_t& ctau, Double_t& ctauErr) {
    TVector3 vtx(decayVrtx->position().x(), decayVrtx->position().y(), 0);
    TVector3 pvtx(bs.position().x(), bs.position().y(), 0);
    TVector3 pperp(kinePart->currentState().globalMomentum().x(), kinePart->currentState().globalMomentum().y(), 0);
    TVector3 vdiff = vtx - pvtx;
    double LxyPV = vdiff.Dot(pperp) / pperp.Mag();
    ctau = LxyPV * kinePart->currentState().mass() / pperp.Perp();

    GlobalError DecayErr = decayVrtx->error();
    GlobalError PrimaryErr = bs.error();
    AlgebraicVector vpperp(3);
    vpperp[0] = pperp.x();
    vpperp[1] = pperp.y();
    vpperp[2] = 0.;
    AlgebraicSymMatrix vXYe = asHepMatrix(DecayErr.matrix()) + asHepMatrix(PrimaryErr.matrix());
    ctauErr = sqrt(vXYe.similarity(vpperp)) * kinePart->currentState().mass() / (pperp.Perp2());
    
    return;
}

// define this as a plug-in
DEFINE_FWK_MODULE(NtupleMaker);
