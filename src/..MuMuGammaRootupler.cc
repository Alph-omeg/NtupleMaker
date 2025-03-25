// -*- C++ -*-
//
// Package:    MuMuGammaRootupler
// Class:      MuMuGammaRootupler
// 
// Description: Dump mumugamma decays
//
// Author:  Zhen Hu
//

// system include files
#include <memory>
//#include "../interface/MultiLepPAT.h"
#include "../interface/VertexReProducer.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
//#include "DataFormats/VertexReco/interface/VertexFwd.h"
//For kinematic fit:
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
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

#include "DataFormats/PatCandidates/interface/Muon.h"
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/Common/interface/View.h>
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH2F.h"
//GEN
//10.10
#include "HeavyFlavorAnalysis/Onia2MuMu//interface/HZZ4LGENAna.h"
//#include "FourMuonAna/Onia/interface/HZZ4LGENAna.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace muon;
using namespace trigger;

std::vector<std::vector<pat::Muon>> muons_previousEvent;
std::vector<pat::Muon> muons_previousEvent_bestYMass;
const ParticleMass muonMass(0.10565837);
float muonSigma = muonMass*1E-6;

//Jinfeng
Float_t temp_Lxy_noMC;
Float_t temp_Lxy_MC;
Float_t temp_CosAlpha_noMC;
Float_t temp_CosAlpha_MC;
Float_t temp_LxyPV_noMC;
Float_t temp_LxyPV_MC;

//Jinfeng 10.3
Float_t temp_LxyPVErr_noMC;
Float_t temp_LxyPVErr_MC;
Float_t temp_LxyPVSig_noMC;
Float_t temp_LxyPVSig_MC;
//The 'Distance' here means the distance between two Jpsi vertexes
Float_t temp_Distance_noMC;
Float_t temp_Distance_MC;
Float_t temp_DistanceErr_noMC;
Float_t temp_DistanceErr_MC;
Float_t temp_DistanceSig_noMC;
Float_t temp_DistanceSig_MC;

//Jinfeng 12.14
bool Dimuon_trigger_matched_1 = false;
bool Dimuon_trigger_matched_2 = false;
bool Rest_Muon_trigger_Matched_1 = false;
bool Rest_Muon_trigger_Matched_2 = false;
bool Trigger_Matched = false;

//Jinfeng 1.16
float DuplicateHLT_Pt_1;
float DuplicateHLT_Pt_2;


//
// class declaration
//

class MuMuGammaRootupler:public edm::EDAnalyzer {
	public:
		explicit MuMuGammaRootupler(const edm::ParameterSet &);
		~MuMuGammaRootupler();

		static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

	private:
		UInt_t getTriggerBits(const edm::Event &);

        //Jinfeng 12.14
        //Modify the input of this function
        bool TriggerMatch(TLorentzVector mu3p4, TLorentzVector mu4p4);
        bool TriggerMatch_restMuons(TLorentzVector mu3p4, TLorentzVector mu4p4);
        bool findTrigger(edm::Handle<edm::TriggerResults>& hltR,
            std::vector < std::string >& triggersToCheck,
            std::vector < std::string >& triggerNameFound);
        void analyzeTrigger(edm::Handle<edm::TriggerResults>& hltR,
            edm::Handle<trigger::TriggerEvent>& hltE,
            const std::string& triggerName);
        bool triggerDecision(edm::Handle<edm::TriggerResults>& hltR, int iTrigger);
		bool   isAncestor(const reco::Candidate *, const reco::Candidate *);
		const  reco::Candidate* GetAncestor(const reco::Candidate *);
		int   tightMuon(edm::View<pat::Muon>::const_iterator rmu, reco::Vertex vertex);
		int   mediumMuon(edm::View<pat::Muon>::const_iterator rmu);
		int   looseMuon(edm::View<pat::Muon>::const_iterator rmu);
		void fillUpsilonVector(RefCountedKinematicTree mumuVertexFitTree, pat::CompositeCandidate dimuonCand, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs);
		void fillUpsilonBestVertex(RefCountedKinematicTree mumuVertexFitTree, pat::CompositeCandidate dimuonCand, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs);
		void fillUpsilonBestMass(RefCountedKinematicTree mumuVertexFitTree, pat::CompositeCandidate dimuonCand, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs);
		void fourMuonFit(pat::CompositeCandidate dimuonCand, RefCountedKinematicTree treeupsilon_part1, edm::Handle< edm::View<pat::Muon> > muons, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV);
		int fourMuonMixFit(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, std::vector<pat::Muon> muons_previous, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV);
		int fourMuonMixFit(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, std::vector<pat::Muon> muons_previous1, std::vector<pat::Muon> muons_previous2, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV);      
		void fourMuonFit_bestYMass(pat::CompositeCandidate dimuonCand, RefCountedKinematicTree treeupsilon_part1, edm::Handle< edm::View<pat::Muon> > muons, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV);
                double GetcTau(RefCountedKinematicVertex& decayVrtx, RefCountedKinematicParticle& kinePart, Vertex& bs, bool MassConstraint);
                double GetcTauErr(RefCountedKinematicVertex& decayVrtx, RefCountedKinematicParticle& kinePart, Vertex& bs);

                //Jinfeng 10.3
                void GetDistance(RefCountedKinematicVertex decayVrtx1, RefCountedKinematicVertex decayVrtx2, bool MassConstraint);
		int  fourMuonMixFit_bestYMass(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, std::vector<pat::Muon> muons_previous, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV);
               //Stand Alone Function for YY
                void YY_fourMuonFit(edm::Handle< edm::View<pat::Muon> > muons,edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV,edm::Handle<reco::VertexCollection>);
                //void setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles, edm::Handle<pat::PackedGenParticle> packedgenParticles);
                void setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles);
		virtual void beginJob();
		virtual void analyze(const edm::Event &, const edm::EventSetup &);
		virtual void endJob();

		virtual void beginRun(edm::Run &, edm::EventSetup const&);
		virtual void endRun(edm::Run const &, edm::EventSetup const &);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);
		virtual void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);
                //GEN
                HZZ4LGENAna genAna;
		// ----------member data ---------------------------
		std::string file_name;
                std::string triggersPassed;
                vector<string> mu1_filtersMatched; // for mu1, all filters it is matched to
                vector<string> mu2_filtersMatched; // for mu2, all filters it is matched to
		edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_Label;
		edm::EDGetTokenT<pat::CompositeCandidateCollection> conversion_Label;
		edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
		edm::EDGetTokenT<reco::BeamSpot> bs_Label;
		edm::EDGetTokenT<edm::View<pat::Muon>> muon_Label;
		edm::EDGetTokenT<edm::TriggerResults> triggerResultsTok_;
                edm::EDGetTokenT<trigger::TriggerEvent>triggerEventTok_;
//                edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
		int  pdgid_;
                vector<string> triggersToApply;
		std::vector<double> OniaMassCuts_;
		bool isMC_;
                bool run_YY_SAF_;
                bool isTriggerMatch_;
		bool OnlyBest_;
		bool OnlyGen_;
		double upsilon_mass_;
                double upsilon_sigma_;
                double upsilon_vtx_cut;
		uint32_t triggerCuts_;
		bool best4muonCand_;

		TH1F* myFourMuM_fit;
		TH1F* myFourMuVtxP_fit;

		TH1F* myDimuonMass_all;
 		TH1F* myY1SFitMass_all;
                TH1F* Total_events;
                int Total_events_analyzed;
	        int Total_events_analyzed_triggered;
	        int Total_events_dimuon;
	        int Total_events_dimuon_trg_matched;
                int Total_events_Fourmuon;
                int nGoodUpsilonCand ;
                int nTriggeredUpsilonCand ;
                int nTriggeredFourMuonCand ;
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
		Int_t  irank;
		UInt_t trigger;
		Int_t numPrimaryVertices;
		Float_t pv_x;
		Float_t pv_y;
		Float_t pv_z;
                Float_t ups1_mass_GenMatched;
                Float_t ups2_mass_GenMatched;
                Float_t ups1_y_GenMatched;
                Float_t ups2_y_GenMatched;
                Float_t ups1_pt_GenMatched;
                Float_t ups2_pt_GenMatched;
                Float_t ups1_mass_GenMatched_ID;
                Float_t ups1_pt_GenMatched_ID;
                Float_t ups1_y_GenMatched_ID;
                Float_t ups2_mass_GenMatched_ID;
                Float_t ups2_pt_GenMatched_ID;
                Float_t ups2_y_GenMatched_ID;
                Float_t ups1_mass_GenMatched_ID_OS_VTX;
                Float_t ups1_pt_GenMatched_ID_OS_VTX;
                Float_t ups1_y_GenMatched_ID_OS_VTX;
                Float_t ups2_mass_GenMatched_ID_OS_VTX;
                Float_t ups2_pt_GenMatched_ID_OS_VTX;
                Float_t ups2_y_GenMatched_ID_OS_VTX;
                Float_t fourmu_mass_allcuts;
                Float_t fourmu_mass_trigger;

                

		std::vector<Int_t> v_mu1Charge;
		std::vector<Int_t> v_mu2Charge;
		std::vector<Float_t> v_mu1_d0;
		std::vector<Float_t> v_mu2_d0;
		std::vector<Float_t> v_mu1_d0err;
		std::vector<Float_t> v_mu2_d0err;
		std::vector<Float_t> v_mu1_dz;
		std::vector<Float_t> v_mu2_dz;
		std::vector<Float_t> v_mu1_dzerr;
		std::vector<Float_t> v_mu2_dzerr;
		std::vector<Float_t> v_mu1_vz;
		std::vector<Float_t> v_mu2_vz;
                std::vector<Float_t> v_mu1_trg_dR;
                std::vector<Float_t> v_mu2_trg_dR;
		std::vector<Float_t> v_mumufit_Mass;
		std::vector<Float_t> v_mumufit_MassErr;
                std::vector<Float_t> v_mumufit_Mass_noMC;
                std::vector<Float_t> v_mumufit_MassErr_noMC;
                std::vector<Float_t> v_mumu_Pt;
                std::vector<Float_t> v_mumu_Rapidity;    
		std::vector<Float_t> v_mumufit_VtxCL;
                std::vector<Float_t> v_mumufit_VtxCL_noMC;
		std::vector<Float_t> v_mumufit_VtxCL2;
                std::vector<Float_t> v_mumufit_cTau_noMC;
                std::vector<Float_t> v_mumufit_cTau_MC;
                std::vector<Float_t> v_mumufit_cTauErr_noMC;
                std::vector<Float_t> v_mumufit_cTauErr_MC;

        //Jinfeng
        std::vector<Float_t> v_mumufit_Lxy_noMC;
        std::vector<Float_t> v_mumufit_Lxy_MC;
        std::vector<Float_t> v_mumufit_CosAlpha_noMC;
        std::vector<Float_t> v_mumufit_CosAlpha_MC;
        std::vector<Float_t> v_mumufit_LxyPV_noMC;
        std::vector<Float_t> v_mumufit_LxyPV_MC;

        //Jinfeng 10.3
        std::vector<Float_t> v_mumufit_LxyPVErr_noMC;
        std::vector<Float_t> v_mumufit_LxyPVErr_MC;
        std::vector<Float_t> v_mumufit_LxyPVSig_noMC;
        std::vector<Float_t> v_mumufit_LxyPVSig_MC;

        std::vector<RefCountedKinematicVertex> Vector_DimuonVertex_noMC;
        std::vector<RefCountedKinematicVertex> Vector_DimuonVertex_MC;

		std::vector<Double_t> v_mumufit_DecayVtxX;
		std::vector<Double_t> v_mumufit_DecayVtxY;
		std::vector<Double_t> v_mumufit_DecayVtxZ;
		std::vector<Double_t> v_mumufit_DecayVtxXE;
		std::vector<Double_t> v_mumufit_DecayVtxYE;
		std::vector<Double_t> v_mumufit_DecayVtxZE;

		TLorentzVector dimuon_p4;
		TLorentzVector mu1_p4;
		TLorentzVector mu2_p4;
                Float_t mu1charge;
                Float_t mu2charge;
		Float_t mu1_d0;
		Float_t mu2_d0;
		Float_t mu1_d0err;
		Float_t mu2_d0err;
		Float_t mu1_dz;
		Float_t mu2_dz;
		Float_t mu1_dzerr;
		Float_t mu2_dzerr;
		Float_t mu1_vz;
		Float_t mu2_vz;
		Float_t mumufit_Mass;
		Float_t mumufit_MassErr;
		Float_t mumufit_VtxCL;
		Float_t mumufit_VtxCL2;
		Double_t mumufit_DecayVtxX;
		Double_t mumufit_DecayVtxY;
		Double_t mumufit_DecayVtxZ;
		Double_t mumufit_DecayVtxXE;
		Double_t mumufit_DecayVtxYE;
		Double_t mumufit_DecayVtxZE;
		TLorentzVector mumufit_p4;

		std::vector<Float_t> fourMuFit_Mass_allComb_mix;
		Float_t fourMuFit_Mass_mix;
		Float_t fourMuFit_MassErr_mix;
		Float_t fourMuFit_VtxX_mix;
		Float_t fourMuFit_VtxY_mix;
		Float_t fourMuFit_VtxZ_mix;
		Float_t fourMuFit_VtxProb_mix;
		Float_t fourMuFit_Chi2_mix;
		TLorentzVector fourMuFit_p4_mix;
		Int_t fourMuFit_ndof_mix;
		Int_t fourMuFit_3plus1_mix;
		Int_t mu3Charge_mix;
		Int_t mu4Charge_mix;
		Float_t mu3_d0_mix;
		Float_t mu4_d0_mix;
		Float_t mu3_d0err_mix;
		Float_t mu4_d0err_mix;
		Float_t mu3_dz_mix;
		Float_t mu4_dz_mix;
		Float_t mu3_dzerr_mix;
		Float_t mu4_dzerr_mix;
		TLorentzVector fourMuFit_mu1p4_mix;
		TLorentzVector fourMuFit_mu2p4_mix;
		TLorentzVector fourMuFit_mu3p4_mix;
		TLorentzVector fourMuFit_mu4p4_mix;
		TLorentzVector mu3_p4_mix;
		TLorentzVector mu4_p4_mix;
		Float_t fourMuFit_Mass_mix3evts;
		Float_t fourMuFit_VtxProb_mix3evts;
		TLorentzVector fourMuFit_p4_mix3evts;

                // -------------------------
                // GEN level information
                // -------------------------
                // Event variables
                int GENfinalState;
                bool passedFiducialSelection;
                std::vector<double> GENmu_pt; std::vector<double> GENmu_eta; std::vector<double> GENmu_phi; std::vector<double> GENmu_mass; 
                std::vector<int> GENmu_id; std::vector<int> GENmu_status; 
                std::vector<int> GENmu_MomId; std::vector<int> GENmu_MomMomId;
                int GENlep_Xindex[4];
                std::vector<int> All_ups_muindex1;
                std::vector<int> All_ups_muindex2;
                int GEN_ups1_mu1_index;
                int GEN_ups1_mu2_index;
                int GEN_ups2_mu1_index;
                int GEN_ups2_mu2_index;
                int RECO_ups1_mu1_index;
                int RECO_ups1_mu2_index;
                int RECO_ups2_mu1_index;
                int RECO_ups2_mu2_index;
                std::vector<double> AllRecoMuons_Pt;
                std::vector<double> AllRecoMuons_Eta;
                std::vector<int> GENups_DaughtersId;
                std::vector<double> GENups_Daughter_mupt;
                std::vector<double> GENups_Daughter_mueta;
                std::vector<double> GENups_Daughter_muphi;
                std::vector<double> GENups_Daughter_mumass;
                std::vector<double> GENups_Daughter_mustatus;
                std::vector<int> GENups_Daughter_muid;
                std::vector<int> GENups_MomId;
                std::vector<double> GENups_pt;
                std::vector<double> GENups_eta;
                std::vector<double> GENups_y;
                std::vector<double> GENups_phi;
                std::vector<double> GENups_mass;
                std::vector<double> GENdimu_mass;
                std::vector<double> GENdimu_pt;
                std::vector<double> GENdimu_eta;
                std::vector<double> GENdimu_y;
                std::vector<double> GENdimu_phi;
                std::vector<double> GENX_mass;
                std::vector<double> GENX_pt;
                std::vector<double> GENX_eta;
                std::vector<double> GENX_y;
                std::vector<double> GENX_phi;
		std::vector<Float_t> fourMuFit_Mass_allComb;
		std::vector<Float_t> fourMuFit_Mass;
		std::vector<Float_t> fourMuFit_MassErr;
		std::vector<Float_t> fourMuFit_Pt;
		std::vector<Float_t> fourMuFit_Eta;
		std::vector<Float_t> fourMuFit_Phi;
		std::vector<Double_t> fourMuFit_VtxX;
		std::vector<Double_t> fourMuFit_VtxY;
		std::vector<Double_t> fourMuFit_VtxZ;
                std::vector<Double_t> fourMuFit_VtxXE;
                std::vector<Double_t> fourMuFit_VtxYE;
                std::vector<Double_t> fourMuFit_VtxZE;
                std::vector<Double_t> fourMuFit_PVX;
                std::vector<Double_t> fourMuFit_PVY;
                std::vector<Double_t> fourMuFit_PVZ;
                std::vector<Double_t> fourMuFit_PVXE;
                std::vector<Double_t> fourMuFit_PVYE;
                std::vector<Double_t> fourMuFit_PVZE;
                std::vector<Double_t> fourMuFit_PV_XYE;
                std::vector<Double_t> fourMuFit_PV_XZE;
                std::vector<Double_t> fourMuFit_PV_YZE;
                std::vector<Double_t> fourMuFit_PV_CL;
		std::vector<Float_t> fourMuFit_VtxProb;

                //Jinfeng_Xdistance
        std::vector<Double_t> fourMuFit_cosAlpha_noMC;
        std::vector<Double_t> fourMuFit_lxy_noMC;
        std::vector<Double_t> fourMuFit_lxyerr_noMC;
        std::vector<Double_t> fourMuFit_lxySig_noMC;
        std::vector<Double_t> fourMuFit_ctauPV_noMC;
        std::vector<Double_t> fourMuFit_ctauErrPV_noMC;
        std::vector<Double_t> fourMuFit_lxyPV_noMC;
        std::vector<Double_t> fourMuFit_cosAlpha;
        std::vector<Double_t> fourMuFit_lxy;
        std::vector<Double_t> fourMuFit_lxyerr;
        std::vector<Double_t> fourMuFit_lxySig;
        std::vector<Double_t> fourMuFit_ctauPV;
        std::vector<Double_t> fourMuFit_ctauErrPV;
        std::vector<Double_t> fourMuFit_lxyPV;

        //Jinfeng 10.3
        std::vector<Float_t> fourMuFit_Distance_noMC;
        std::vector<Float_t> fourMuFit_Distance_MC;
        std::vector<Float_t> fourMuFit_DistanceErr_noMC;
        std::vector<Float_t> fourMuFit_DistanceErr_MC;
        std::vector<Float_t> fourMuFit_DistanceSig_noMC;
        std::vector<Float_t> fourMuFit_DistanceSig_MC;
      
                std::vector<Float_t> fourMuFit_Double_ups_VtxProb;
                std::vector<Float_t> fourMuFit_ups1_VtxProb;
                std::vector<Float_t> fourMuFit_ups2_VtxProb;
                std::vector<Float_t> fourMuFit_ups1_VtxProb_noMC;
                std::vector<Float_t> fourMuFit_ups2_VtxProb_noMC;
                std::vector<Float_t> fourMuFit_ups1_mass;
                std::vector<Float_t> fourMuFit_ups2_mass;
                std::vector<Float_t> fourMuFit_ups1_cTau_noMC;
                std::vector<Float_t> fourMuFit_ups2_cTau_noMC;
                std::vector<Float_t> fourMuFit_ups1_cTau_MC;
                std::vector<Float_t> fourMuFit_ups2_cTau_MC;
                std::vector<Float_t> fourMuFit_ups1_cTauErr_noMC;
                std::vector<Float_t> fourMuFit_ups2_cTauErr_noMC;
                std::vector<Float_t> fourMuFit_ups1_cTauErr_MC;
                std::vector<Float_t> fourMuFit_ups2_cTauErr_MC;

                //Jinfeng
                std::vector<Float_t> fourMuFit_ups1_Lxy_noMC;
                std::vector<Float_t> fourMuFit_ups2_Lxy_noMC;
                std::vector<Float_t> fourMuFit_ups1_Lxy_MC;
                std::vector<Float_t> fourMuFit_ups2_Lxy_MC;
                std::vector<Float_t> fourMuFit_ups1_CosAlpha_noMC;
                std::vector<Float_t> fourMuFit_ups2_CosAlpha_noMC;
                std::vector<Float_t> fourMuFit_ups1_CosAlpha_MC;
                std::vector<Float_t> fourMuFit_ups2_CosAlpha_MC;

                //Jinfeng 10.3
                std::vector<Float_t> fourMuFit_ups1_LxyPVErr_noMC;
                std::vector<Float_t> fourMuFit_ups2_LxyPVErr_noMC;
                std::vector<Float_t> fourMuFit_ups1_LxyPVErr_MC;
                std::vector<Float_t> fourMuFit_ups2_LxyPVErr_MC;
                std::vector<Float_t> fourMuFit_ups1_LxyPVSig_noMC;
                std::vector<Float_t> fourMuFit_ups2_LxyPVSig_noMC;
                std::vector<Float_t> fourMuFit_ups1_LxyPVSig_MC;
                std::vector<Float_t> fourMuFit_ups2_LxyPVSig_MC;

                std::vector<Float_t> fourMuFit_ups1_LxyPV_noMC;
                std::vector<Float_t> fourMuFit_ups2_LxyPV_noMC;
                std::vector<Float_t> fourMuFit_ups1_LxyPV_MC;
                std::vector<Float_t> fourMuFit_ups2_LxyPV_MC;

                std::vector<Float_t> fourMuFit_ups1_pt;
                std::vector<Float_t> fourMuFit_ups2_pt;
                std::vector<Float_t> fourMuFit_ups1_rapidity;
                std::vector<Float_t> fourMuFit_ups2_rapidity;
                std::vector<Float_t> fourMuFit_ups1_massError;
                std::vector<Float_t> fourMuFit_ups2_massError;
                std::vector<Float_t> fourMuFit_wrong_ups1_mass;
                std::vector<Float_t> fourMuFit_wrong_ups2_mass;
                std::vector<Float_t> fourMuFit_wrong_ups1_massError;
                std::vector<Float_t> fourMuFit_wrong_ups2_massError;
                std::vector<Double_t> fourMuFit_ups1_VtxX;
                std::vector<Double_t> fourMuFit_ups1_VtxY;
                std::vector<Double_t> fourMuFit_ups1_VtxZ;
                std::vector<Double_t> fourMuFit_ups1_VtxXE;
                std::vector<Double_t> fourMuFit_ups1_VtxYE;
                std::vector<Double_t> fourMuFit_ups1_VtxZE;               
                std::vector<Double_t> fourMuFit_ups2_VtxX;
                std::vector<Double_t> fourMuFit_ups2_VtxY;
                std::vector<Double_t> fourMuFit_ups2_VtxZ;
                std::vector<Double_t> fourMuFit_ups2_VtxXE;
                std::vector<Double_t> fourMuFit_ups2_VtxYE;
                std::vector<Double_t> fourMuFit_ups2_VtxZE;                           
		std::vector<Float_t> fourMuFit_Chi2;
		std::vector<Int_t> fourMuFit_ndof;
                std::vector<Float_t> fourMuFit_Double_ups_Chi2;
                std::vector<Int_t> fourMuFit_Double_ups_ndof;
		std::vector<Float_t> fourMuFit_mu1Pt;
		std::vector<Float_t> fourMuFit_mu1Eta;
		std::vector<Float_t> fourMuFit_mu1Phi;
		std::vector<Float_t> fourMuFit_mu1E;
                std::vector<Int_t> fourMuFit_mu1frompv;
		std::vector<Float_t> fourMuFit_mu2Pt;
		std::vector<Float_t> fourMuFit_mu2Eta;
		std::vector<Float_t> fourMuFit_mu2Phi;
		std::vector<Float_t> fourMuFit_mu2E;
                std::vector<Int_t> fourMuFit_mu2frompv;
		std::vector<Float_t> fourMuFit_mu3Pt;
		std::vector<Float_t> fourMuFit_mu3Eta;
		std::vector<Float_t> fourMuFit_mu3Phi;
		std::vector<Float_t> fourMuFit_mu3E;
                std::vector<Int_t> fourMuFit_mu3frompv;
		std::vector<Float_t> fourMuFit_mu4Pt;
		std::vector<Float_t> fourMuFit_mu4Eta;
		std::vector<Float_t> fourMuFit_mu4Phi;
		std::vector<Float_t> fourMuFit_mu4E;
                std::vector<Int_t> fourMuFit_mu4frompv;
                std::vector<Float_t> fourMuFit_mu3_trg_dR;
                std::vector<Float_t> fourMuFit_mu4_trg_dR;
                std::vector<Float_t> mu1_Pt;
                std::vector<Float_t> mu1_Eta;
                std::vector<Float_t> mu1_Phi;
                std::vector<Float_t> mu1_E;
                std::vector<int> mu1_genindex;
                std::vector<Float_t> mu2_Pt;
                std::vector<Float_t> mu2_Eta;
                std::vector<Float_t> mu2_Phi;
                std::vector<Float_t> mu2_E;
                std::vector<int> mu2_genindex;
		std::vector<Float_t> mu3_Pt;
		std::vector<Float_t> mu3_Eta;
		std::vector<Float_t> mu3_Phi;
		std::vector<Float_t> mu3_E;
		std::vector<Float_t> mu4_Pt;
                //std::vector<int> mu3_genindex;
		std::vector<Float_t> mu4_Eta;
		std::vector<Float_t> mu4_Phi;
		std::vector<Float_t> mu4_E;
                //std::vector<int> mu4_genindex;
                std::vector<Int_t> mu1Charge;
                std::vector<Int_t> mu2Charge;
		std::vector<Int_t> mu3Charge;
		std::vector<Int_t> mu4Charge;
		std::vector<Float_t> mu3_d0;
		std::vector<Float_t> mu4_d0;
		std::vector<Float_t> mu3_d0err;
		std::vector<Float_t> mu4_d0err;
		std::vector<Float_t> mu3_dz;
		std::vector<Float_t> mu4_dz;
		std::vector<Float_t> mu3_dzerr;
		std::vector<Float_t> mu4_dzerr;
                std::vector<Int_t> mu1_Soft;
                std::vector<Int_t> mu2_Soft;
                std::vector<Int_t> mu3_Soft;
                std::vector<Int_t> mu4_Soft;
		std::vector<Int_t> mu1_Tight;
		std::vector<Int_t> mu2_Tight;
		std::vector<Int_t> mu3_Tight;
		std::vector<Int_t> mu4_Tight;
		std::vector<Int_t> mu1_Medium;
		std::vector<Int_t> mu2_Medium;
		std::vector<Int_t> mu3_Medium;
		std::vector<Int_t> mu4_Medium;
		std::vector<Int_t> mu1_Loose;
		std::vector<Int_t> mu2_Loose;
		std::vector<Int_t> mu3_Loose;
		std::vector<Int_t> mu4_Loose;
		std::vector<Int_t> mu1_pdgID;
		std::vector<Int_t> mu2_pdgID;
		std::vector<Int_t> mu3_pdgID;
		std::vector<Int_t> mu4_pdgID;
                std::vector<Int_t> fourMuFit_mu12overlap;
                std::vector<Int_t> fourMuFit_mu13overlap;
                std::vector<Int_t> fourMuFit_mu14overlap;
                std::vector<Int_t> fourMuFit_mu23overlap;
                std::vector<Int_t> fourMuFit_mu24overlap;
                std::vector<Int_t> fourMuFit_mu34overlap;
                std::vector<Int_t> fourMuFit_mu12SharedSeg;
                std::vector<Int_t> fourMuFit_mu13SharedSeg;
                std::vector<Int_t> fourMuFit_mu14SharedSeg;
                std::vector<Int_t> fourMuFit_mu23SharedSeg;
                std::vector<Int_t> fourMuFit_mu24SharedSeg;
                std::vector<Int_t> fourMuFit_mu34SharedSeg;
		/*
			std::vector<Float_t> fourMuFit_Mass_allComb;
			Float_t fourMuFit_Mass;
			Float_t fourMuFit_MassErr;
			Float_t fourMuFit_VtxX;
			Float_t fourMuFit_VtxY;
			Float_t fourMuFit_VtxZ;
			Float_t fourMuFit_VtxProb;
			Float_t fourMuFit_Chi2;
			Int_t fourMuFit_ndof;
			TLorentzVector fourMuFit_p4;
			Int_t mu3Charge;
			Int_t mu4Charge;
			Float_t mu3_d0;
			Float_t mu4_d0;
			Float_t mu3_d0err;
			Float_t mu4_d0err;
			Float_t mu3_dz;
			Float_t mu4_dz;
			Float_t mu3_dzerr;
			Float_t mu4_dzerr;
			TLorentzVector fourMuFit_mu1p4;
			TLorentzVector fourMuFit_mu2p4;
			TLorentzVector fourMuFit_mu3p4;
			TLorentzVector fourMuFit_mu4p4;
			TLorentzVector mu3_p4;
			TLorentzVector mu4_p4;
			Int_t mu1_Tight;
			Int_t mu2_Tight;
			Int_t mu3_Tight;
			Int_t mu4_Tight;
			Int_t mu3_pdgID;
			Int_t mu4_pdgID;
			std::vector<pat::CompositeCandidate> upsilonMuons;
			std::vector<pat::Muon> theRestMuons;
			*/

		TLorentzVector dimuon_p4_bestYMass;
		TLorentzVector mu1_p4_bestYMass;
		TLorentzVector mu2_p4_bestYMass;
		Int_t mu1Charge_bestYMass;
		Int_t mu2Charge_bestYMass;
		Float_t mu1_d0_bestYMass;
		Float_t mu2_d0_bestYMass;
		Float_t mu1_d0err_bestYMass;
		Float_t mu2_d0err_bestYMass;
		Float_t mu1_dz_bestYMass;
		Float_t mu2_dz_bestYMass;
		Float_t mu1_dzerr_bestYMass;
		Float_t mu2_dzerr_bestYMass;
		Float_t mumufit_Mass_bestYMass;
		Float_t mumufit_MassErr_bestYMass;
		Float_t mumufit_VtxCL_bestYMass;
		Float_t mumufit_VtxCL2_bestYMass;
		Float_t mumufit_DecayVtxX_bestYMass;
		Float_t mumufit_DecayVtxY_bestYMass;
		Float_t mumufit_DecayVtxZ_bestYMass;
		Float_t mumufit_DecayVtxXE_bestYMass;
		Float_t mumufit_DecayVtxYE_bestYMass;
		Float_t mumufit_DecayVtxZE_bestYMass;
		TLorentzVector mumufit_p4_bestYMass;
		Int_t bestVertex_and_bestYMass;

		std::vector<Float_t> fourMuFit_Mass_allComb_mix_bestYMass;
		Float_t fourMuFit_Mass_mix_bestYMass;
		Float_t fourMuFit_MassErr_mix_bestYMass;
		Float_t fourMuFit_VtxX_mix_bestYMass;
		Float_t fourMuFit_VtxY_mix_bestYMass;
		Float_t fourMuFit_VtxZ_mix_bestYMass;
		Float_t fourMuFit_VtxProb_mix_bestYMass;
		Float_t fourMuFit_Chi2_mix_bestYMass;
		Int_t fourMuFit_ndof_mix_bestYMass;
		Int_t fourMuFit_3plus1_mix_bestYMass;
		TLorentzVector fourMuFit_p4_mix_bestYMass;
		Int_t mu3Charge_mix_bestYMass;
		Int_t mu4Charge_mix_bestYMass;
		Float_t mu3_d0_mix_bestYMass;
		Float_t mu4_d0_mix_bestYMass;
		Float_t mu3_d0err_mix_bestYMass;
		Float_t mu4_d0err_mix_bestYMass;
		Float_t mu3_dz_mix_bestYMass;
		Float_t mu4_dz_mix_bestYMass;
		Float_t mu3_dzerr_mix_bestYMass;
		Float_t mu4_dzerr_mix_bestYMass;
		TLorentzVector fourMuFit_mu1p4_mix_bestYMass;
		TLorentzVector fourMuFit_mu2p4_mix_bestYMass;
		TLorentzVector fourMuFit_mu3p4_mix_bestYMass;
		TLorentzVector fourMuFit_mu4p4_mix_bestYMass;
		TLorentzVector mu3_p4_mix_bestYMass;
		TLorentzVector mu4_p4_mix_bestYMass;

		std::vector<Float_t> fourMuFit_Mass_allComb_bestYMass;
		Float_t fourMuFit_Mass_bestYMass;
		Float_t fourMuFit_MassErr_bestYMass;
		Float_t fourMuFit_VtxX_bestYMass;
		Float_t fourMuFit_VtxY_bestYMass;
		Float_t fourMuFit_VtxZ_bestYMass;
		Float_t fourMuFit_VtxProb_bestYMass;
		Float_t fourMuFit_Chi2_bestYMass;
		Int_t fourMuFit_ndof_bestYMass;
		TLorentzVector fourMuFit_p4_bestYMass;
		Int_t mu3Charge_bestYMass;
		Int_t mu4Charge_bestYMass;
		Float_t mu3_d0_bestYMass;
		Float_t mu4_d0_bestYMass;
		Float_t mu3_d0err_bestYMass;
		Float_t mu4_d0err_bestYMass;
		Float_t mu3_dz_bestYMass;
		Float_t mu4_dz_bestYMass;
		Float_t mu3_dzerr_bestYMass;
		Float_t mu4_dzerr_bestYMass;
		TLorentzVector fourMuFit_mu1p4_bestYMass;
		TLorentzVector fourMuFit_mu2p4_bestYMass;
		TLorentzVector fourMuFit_mu3p4_bestYMass;
		TLorentzVector fourMuFit_mu4p4_bestYMass;
		TLorentzVector mu3_p4_bestYMass;
		TLorentzVector mu4_p4_bestYMass;
		Int_t mu1_Tight_bestYMass;
		Int_t mu2_Tight_bestYMass;
		Int_t mu3_Tight_bestYMass;
		Int_t mu4_Tight_bestYMass;
		Int_t mu3_pdgID_bestYMass;
		Int_t mu4_pdgID_bestYMass;
                Int_t mu1_Medium_bestYMass; 
                Int_t mu2_Medium_bestYMass;
                Int_t mu3_Medium_bestYMass;
                Int_t mu4_Medium_bestYMass;  
                Int_t mu1_Loose_bestYMass;
                Int_t mu2_Loose_bestYMass;
                Int_t mu3_Loose_bestYMass;
                Int_t mu4_Loose_bestYMass;
                std::string rootFileName;
                edm::InputTag triggerEventTag_;   
                std::string hltName_;
                std::string   triggerName_;
                std::string hlTriggerSummaryAOD_; 
                edm::TriggerNames triggerNames;
                bool verbose;
                float trg_Match_dR_cut, trg_Match_dP_cut,trg_Match_dP_ratio_cut; 
                std::vector<std::string> triggerList;
                bool checkTrigger;
                int runNumber;
  		std::vector < reco::MuonCollection::const_iterator > allL1TrigMuons;
  		std::vector < reco::MuonCollection::const_iterator > allL2TrigMuons;
  		std::vector < reco::MuonCollection::const_iterator > allL3TrigMuons;
   		std::vector<TLorentzVector> allTrigMuons;
                std::vector<TLorentzVector> allRestTrigMuons;
  		std::vector < GlobalVector > allMuL1TriggerVectors;
  		std::vector < GlobalVector > allMuL2TriggerVectors;
  		std::vector < GlobalVector > allMuL3TriggerVectors_lowEff;
  		std::vector < GlobalVector > allMuL3TriggerVectors_highEff;
  		std::vector < GlobalVector > allMuHLTTriggerVectors;
                std::vector < GlobalVector > allUpsilonMuHLTTriggerVectors;
                int lastTriggerModule;
                HLTConfigProvider hltConfig_;
		TTree *onia_tree;
		//TTree *gen_tree;

		Int_t mother_pdgId;
		Int_t dimuon_pdgId;
		TLorentzVector gen_dimuon_p4;
		TLorentzVector gen_mu1_p4;
		TLorentzVector gen_mu2_p4;
                edm::EDGetTokenT<reco::GenParticleCollection> prunedgenParticlesSrc_;
                //edm::EDGetTokenT<pat::PackedGenParticle> packedgenParticlesSrc_;
               // Counters
};

//
// constructors and destructor
//
MuMuGammaRootupler::MuMuGammaRootupler(const edm::ParameterSet & iConfig):
	dimuon_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("dimuons"))),
	conversion_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("conversions"))),
	primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
	bs_Label(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("offlineBeamSpot"))),
	muon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter< edm::InputTag>("muons"))),
	triggerResultsTok_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
        triggerEventTok_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("TriggerSummaryAOD"))),
	pdgid_(iConfig.getParameter<uint32_t>("onia_pdgid")),
	OniaMassCuts_(iConfig.getParameter<std::vector<double>>("onia_mass_cuts")),
	isMC_(iConfig.getParameter<bool>("isMC")),
        run_YY_SAF_(iConfig.getParameter<bool>("run_YY_SAF")),
        isTriggerMatch_(iConfig.getParameter<bool>("TriggerMatch")),
	OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
	OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
	upsilon_mass_(iConfig.getParameter<double>("upsilon_mass")),
        upsilon_sigma_(iConfig.getParameter<double>("upsilon_sigma")),
        upsilon_vtx_cut(iConfig.getParameter<double>("upsilon_vertex_CL_cut")),
	triggerCuts_(iConfig.getParameter<uint32_t>("triggerCuts")),
	best4muonCand_(iConfig.getParameter<bool>("best4muonCand")),
        verbose(iConfig.getUntrackedParameter<bool>("VERBOSE",false)),
        trg_Match_dR_cut(iConfig.getParameter<double>("TRG_Match_DR")),
        trg_Match_dP_cut(iConfig.getParameter<double>("TRG_Match_DP")),
        trg_Match_dP_ratio_cut(iConfig.getParameter<double>("TRG_Match_DP_ratio")),
        triggerList(iConfig.getUntrackedParameter<std::vector<std::string>>("triggerList"))
{
	edm::Service < TFileService > fs;
	onia_tree = fs->make < TTree > ("oniaTree", "Tree of MuMuGamma");
	//gen_tree = fs->make < TTree > ("genTree", "Tree of genCand");
        Total_events_analyzed = 0;
        nevent_soft4mu = 0;
        nevent_4mu_pt=0;
        nevent_4mu_eta=0;
	Total_events_analyzed_triggered = 0 ;
	Total_events_dimuon= 0;
	Total_events_dimuon_trg_matched = 0;
	Total_events_Fourmuon = 0;

        runNumber = -99;
        checkTrigger =true;
        hltName_ = "HLT";
        triggerName_ = "@"; // "@" means: analyze all triggers in config
	if (!OnlyGen_) {
		onia_tree->Branch("run",     &run,     "run/l");
		onia_tree->Branch("lumi",     &lumi,     "lumi/l");
		onia_tree->Branch("event",   &event,   "event/l");
		onia_tree->Branch("irank",   &irank,   "irank/I");
		onia_tree->Branch("trigger", &trigger, "trigger/I");
		onia_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
		onia_tree->Branch("pv_x",     &pv_x,     "pv_x/F");
		onia_tree->Branch("pv_y",     &pv_y,     "pv_y/F");
		onia_tree->Branch("pv_z",     &pv_z,     "pv_z/F");
                onia_tree->Branch("ups1_mass_GenMatched", &ups1_mass_GenMatched,"ups1_mass_GenMatched/F");
                onia_tree->Branch("ups2_mass_GenMatched", &ups2_mass_GenMatched,"ups2_mass_GenMatched/F");
                onia_tree->Branch("ups1_y_GenMatched", &ups1_y_GenMatched,"ups1_y_GenMatched/F");
                onia_tree->Branch("ups2_y_GenMatched", &ups2_y_GenMatched,"ups2_y_GenMatched/F");
                onia_tree->Branch("ups1_pt_GenMatched", &ups1_pt_GenMatched,"ups1_pt_GenMatched/F");
                onia_tree->Branch("ups2_pt_GenMatched", &ups2_pt_GenMatched,"ups2_pt_GenMatched/F");
                onia_tree->Branch("ups1_mass_GenMatched_ID", &ups1_mass_GenMatched_ID,"ups1_mass_GenMatched_ID/F");
                onia_tree->Branch("ups1_pt_GenMatched_ID", &ups1_pt_GenMatched_ID,"ups1_pt_GenMatched_ID/F");
                onia_tree->Branch("ups1_y_GenMatched_ID", &ups1_y_GenMatched_ID,"ups1_y_GenMatched_ID/F");
                onia_tree->Branch("ups2_mass_GenMatched_ID", &ups2_mass_GenMatched_ID,"ups2_mass_GenMatched_ID/F");
                onia_tree->Branch("ups2_pt_GenMatched_ID", &ups2_pt_GenMatched_ID,"ups2_pt_GenMatched_ID/F");
                onia_tree->Branch("ups2_y_GenMatched_ID", &ups2_y_GenMatched_ID,"ups2_y_GenMatched_ID/F");
                onia_tree->Branch("ups1_mass_GenMatched_ID_OS_VTX", &ups1_mass_GenMatched_ID_OS_VTX,"ups1_mass_GenMatched_ID_OS_VTX/F");
                onia_tree->Branch("ups1_pt_GenMatched_ID_OS_VTX", &ups1_pt_GenMatched_ID_OS_VTX,"ups1_pt_GenMatched_ID_OS_VTX/F");
                onia_tree->Branch("ups1_y_GenMatched_ID_OS_VTX", &ups1_y_GenMatched_ID_OS_VTX,"ups1_y_GenMatched_ID_OS_VTX/F");
                onia_tree->Branch("ups2_mass_GenMatched_ID_OS_VTX", &ups2_mass_GenMatched_ID_OS_VTX,"ups2_mass_GenMatched_ID_OS_VTX/F");
                onia_tree->Branch("ups2_pt_GenMatched_ID_OS_VTX", &ups2_pt_GenMatched_ID_OS_VTX,"ups2_pt_GenMatched_ID_OS_VTX/F");
                onia_tree->Branch("ups2_y_GenMatched_ID_OS_VTX", &ups2_y_GenMatched_ID_OS_VTX,"ups2_y_GenMatched_ID_OS_VTX/F");
                onia_tree->Branch("fourmu_mass_allcuts", &fourmu_mass_allcuts,"fourmu_mass_allcuts/F");
                onia_tree->Branch("fourmu_mass_trigger", &fourmu_mass_trigger,"fourmu_mass_trigger/F");
                onia_tree->Branch("All_ups_muindex1",&All_ups_muindex1);
                onia_tree->Branch("All_ups_muindex2",&All_ups_muindex2);
                onia_tree->Branch("GEN_ups1_mu1_index",&GEN_ups1_mu1_index,"GEN_ups1_mu1_index/I");
                onia_tree->Branch("GEN_ups1_mu2_index",&GEN_ups1_mu2_index,"GEN_ups1_mu2_index/I");
                onia_tree->Branch("GEN_ups2_mu1_index",&GEN_ups2_mu1_index,"GEN_ups2_mu1_index/I");
                onia_tree->Branch("GEN_ups2_mu2_index",&GEN_ups2_mu2_index,"GEN_ups2_mu2_index/I");
                onia_tree->Branch("RECO_ups1_mu1_index",&RECO_ups1_mu1_index,"RECO_ups1_mu1_index/I");
                onia_tree->Branch("RECO_ups1_mu2_index",&RECO_ups1_mu2_index,"RECO_ups1_mu2_index/I");
                onia_tree->Branch("RECO_ups2_mu1_index",&RECO_ups2_mu1_index,"RECO_ups2_mu1_index/I");
                onia_tree->Branch("RECO_ups2_mu2_index",&RECO_ups2_mu2_index,"RECO_ups2_mu2_index/I");
		onia_tree->Branch("v_mu1Charge",   &v_mu1Charge);
		onia_tree->Branch("v_mu2Charge",   &v_mu2Charge);
		onia_tree->Branch("v_mu1_d0",   &v_mu1_d0);
		onia_tree->Branch("v_mu1_d0err",   &v_mu1_d0err);
		onia_tree->Branch("v_mu2_d0",   &v_mu2_d0);
		onia_tree->Branch("v_mu2_d0err",   &v_mu2_d0err);
		onia_tree->Branch("v_mu1_dz",   &v_mu1_dz);
		onia_tree->Branch("v_mu1_dzerr",   &v_mu1_dzerr);
		onia_tree->Branch("v_mu2_dz",   &v_mu2_dz);
		onia_tree->Branch("v_mu2_dzerr",   &v_mu2_dzerr);
		onia_tree->Branch("v_mu1_vz",   &v_mu1_vz);
		onia_tree->Branch("v_mu2_vz",   &v_mu2_vz);
                onia_tree->Branch("mu1_trg_dR",   &v_mu1_trg_dR);
                onia_tree->Branch("mu2_trg_dR",   &v_mu2_trg_dR);
		onia_tree->Branch("v_mumufit_Mass",&v_mumufit_Mass);
		onia_tree->Branch("v_mumufit_MassErr",&v_mumufit_MassErr);
                onia_tree->Branch("v_mumufit_Mass_noMC",&v_mumufit_Mass_noMC);
                onia_tree->Branch("v_mumufit_MassErr_noMC",&v_mumufit_MassErr_noMC);
                onia_tree->Branch("v_mumu_Pt",&v_mumu_Pt);
                onia_tree->Branch("v_mumu_Rapidity",&v_mumu_Rapidity);
		onia_tree->Branch("v_mumufit_VtxCL",&v_mumufit_VtxCL);
                onia_tree->Branch("v_mumufit_VtxCL_noMC",&v_mumufit_VtxCL_noMC);
		onia_tree->Branch("v_mumufit_VtxCL2",&v_mumufit_VtxCL2);
                onia_tree->Branch("v_mumufit_cTau_noMC",&v_mumufit_cTau_noMC);
                onia_tree->Branch("v_mumufit_cTau_MC",&v_mumufit_cTau_MC);
                onia_tree->Branch("v_mumufit_cTauErr_noMC",&v_mumufit_cTauErr_noMC);
                onia_tree->Branch("v_mumufit_cTauErr_MC",&v_mumufit_cTauErr_MC);

                //Jinfeng
        onia_tree->Branch("v_mumufit_Lxy_noMC", &v_mumufit_Lxy_noMC);
        onia_tree->Branch("v_mumufit_Lxy_MC", &v_mumufit_Lxy_MC);
        onia_tree->Branch("v_mumufit_CosAlpha_noMC", &v_mumufit_CosAlpha_noMC);
        onia_tree->Branch("v_mumufit_CosAlpha_MC", &v_mumufit_CosAlpha_MC);
        onia_tree->Branch("v_mumufit_LxyPV_noMC", &v_mumufit_Lxy_noMC);
        onia_tree->Branch("v_mumufit_LxyPV_MC", &v_mumufit_Lxy_MC);

		onia_tree->Branch("v_mumufit_DecayVtxX",&v_mumufit_DecayVtxX);
		onia_tree->Branch("v_mumufit_DecayVtxY",&v_mumufit_DecayVtxY);
		onia_tree->Branch("v_mumufit_DecayVtxZ",&v_mumufit_DecayVtxZ);
		onia_tree->Branch("v_mumufit_DecayVtxXE",&v_mumufit_DecayVtxXE);
		onia_tree->Branch("v_mumufit_DecayVtxYE",&v_mumufit_DecayVtxYE);
		onia_tree->Branch("v_mumufit_DecayVtxZE",&v_mumufit_DecayVtxZE);
                onia_tree->Branch("mu1_filtersMatched", & mu1_filtersMatched);
                onia_tree->Branch("mu2_filtersMatched", & mu2_filtersMatched);
		onia_tree->Branch("mu1_p4",  "TLorentzVector", &mu1_p4);
		onia_tree->Branch("mu2_p4",  "TLorentzVector", &mu2_p4);
		onia_tree->Branch("mu1charge",   &mu1charge,    "mu1charge/I");
		onia_tree->Branch("mu2charge",   &mu2charge,    "mu2charge/I");
		onia_tree->Branch("mu1_d0",   &mu1_d0,    "mu1_d0/F");
		onia_tree->Branch("mu1_d0err",   &mu1_d0err,    "mu1_d0err/F");
		onia_tree->Branch("mu2_d0",   &mu2_d0,    "mu2_d0/F");
		onia_tree->Branch("mu2_d0err",   &mu2_d0err,    "mu2_d0err/F");
		onia_tree->Branch("mu1_dz",   &mu1_dz,    "mu1_dz/F");
		onia_tree->Branch("mu1_dzerr",   &mu1_dzerr,    "mu1_dzerr/F");
		onia_tree->Branch("mu2_dz",   &mu2_dz,    "mu2_dz/F");
		onia_tree->Branch("mu2_dzerr",   &mu2_dzerr,    "mu2_dzerr/F");
		onia_tree->Branch("mu1_vz",   &mu1_vz,    "mu1_vz/F");
		onia_tree->Branch("mu2_vz",   &mu2_vz,    "mu2_vz/F");
		onia_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
		onia_tree->Branch("mumufit_Mass",&mumufit_Mass,"mumufit_Mass/F");
		onia_tree->Branch("mumufit_MassErr",&mumufit_MassErr,"mumufit_MassErr/F");
		onia_tree->Branch("mumufit_VtxCL",&mumufit_VtxCL,"mumufit_VtxCL/F");
		onia_tree->Branch("mumufit_VtxCL2",&mumufit_VtxCL2,"mumufit_VtxCL2/F");
		onia_tree->Branch("mumufit_DecayVtxX",&mumufit_DecayVtxX,"mumufit_DecayVtxX/D");
		onia_tree->Branch("mumufit_DecayVtxY",&mumufit_DecayVtxY,"mumufit_DecayVtxY/D");
		onia_tree->Branch("mumufit_DecayVtxZ",&mumufit_DecayVtxZ,"mumufit_DecayVtxZ/D");
		onia_tree->Branch("mumufit_DecayVtxXE",&mumufit_DecayVtxXE,"mumufit_DecayVtxXE/D");
		onia_tree->Branch("mumufit_DecayVtxYE",&mumufit_DecayVtxYE,"mumufit_DecayVtxYE/D");
		onia_tree->Branch("mumufit_DecayVtxZE",&mumufit_DecayVtxZE,"mumufit_DecayVtxZE/D");
		onia_tree->Branch("mumufit_p4",  "TLorentzVector", &mumufit_p4);

		onia_tree->Branch("fourMuFit_Mass_allComb_mix",&fourMuFit_Mass_allComb_mix);
		onia_tree->Branch("fourMuFit_Mass_mix",&fourMuFit_Mass_mix,"fourMuFit_Mass_mix/F");
		onia_tree->Branch("fourMuFit_MassErr_mix",&fourMuFit_MassErr_mix,"fourMuFit_MassErr_mix/F");
		onia_tree->Branch("fourMuFit_VtxX_mix",&fourMuFit_VtxX_mix,"fourMuFit_VtxX_mix/F");
		onia_tree->Branch("fourMuFit_VtxY_mix",&fourMuFit_VtxY_mix,"fourMuFit_VtxY_mix/F");
		onia_tree->Branch("fourMuFit_VtxZ_mix",&fourMuFit_VtxZ_mix,"fourMuFit_VtxZ_mix/F");
		onia_tree->Branch("fourMuFit_VtxProb_mix",&fourMuFit_VtxProb_mix,"fourMuFit_VtxProb_mix/F");
		onia_tree->Branch("fourMuFit_Chi2_mix",&fourMuFit_Chi2_mix,"fourMuFit_Chi2_mix/F");
		onia_tree->Branch("fourMuFit_ndof_mix",&fourMuFit_ndof_mix,"fourMuFit_ndof_mix/I");
		onia_tree->Branch("fourMuFit_3plus1_mix",&fourMuFit_3plus1_mix,"fourMuFit_3plus1_mix/I");
		onia_tree->Branch("fourMuFit_p4_mix",  "TLorentzVector", &fourMuFit_p4_mix);
		onia_tree->Branch("fourMuFit_mu1p4_mix",  "TLorentzVector", &fourMuFit_mu1p4_mix);
		onia_tree->Branch("fourMuFit_mu2p4_mix",  "TLorentzVector", &fourMuFit_mu2p4_mix);
		onia_tree->Branch("fourMuFit_mu3p4_mix",  "TLorentzVector", &fourMuFit_mu3p4_mix);
		onia_tree->Branch("fourMuFit_mu4p4_mix",  "TLorentzVector", &fourMuFit_mu4p4_mix);
		onia_tree->Branch("mu3Charge_mix",   &mu3Charge_mix,    "mu3Charge_mix/I");
		onia_tree->Branch("mu4Charge_mix",   &mu4Charge_mix,    "mu4Charge_mix/I");
		onia_tree->Branch("mu3_p4_mix",  "TLorentzVector", &mu3_p4_mix);
		onia_tree->Branch("mu4_p4_mix",  "TLorentzVector", &mu4_p4_mix);
		onia_tree->Branch("mu3_d0_mix",   &mu3_d0_mix,    "mu3_d0_mix/F");
		onia_tree->Branch("mu3_d0err_mix",   &mu3_d0err_mix,    "mu3_d0err_mix/F");
		onia_tree->Branch("mu4_d0_mix",   &mu4_d0_mix,    "mu4_d0_mix/F");
		onia_tree->Branch("mu4_d0err_mix",   &mu4_d0err_mix,    "mu4_d0err_mix/F");
		onia_tree->Branch("mu3_dz_mix",   &mu3_dz_mix,    "mu3_dz_mix/F");
		onia_tree->Branch("mu3_dzerr_mix",   &mu3_dzerr_mix,    "mu3_dzerr_mix/F");
		onia_tree->Branch("mu4_dz_mix",   &mu4_dz_mix,    "mu4_dz_mix/F");
		onia_tree->Branch("mu4_dzerr_mix",   &mu4_dzerr_mix,    "mu4_dzerr_mix/F");
		onia_tree->Branch("fourMuFit_Mass_mix3evts",&fourMuFit_Mass_mix3evts,"fourMuFit_Mass_mix3evts/F");
		onia_tree->Branch("fourMuFit_VtxProb_mix3evts",&fourMuFit_VtxProb_mix3evts,"fourMuFit_VtxProb_mix3evts/F");
		onia_tree->Branch("fourMuFit_p4_mix3evts",  "TLorentzVector", &fourMuFit_p4_mix3evts);
                onia_tree->Branch("GENfinalState",&GENfinalState,"GENfinalState/I");
                onia_tree->Branch("passedFiducialSelection",&passedFiducialSelection,"passedFiducialSelection/O");
                onia_tree->Branch("GENmu_pt",&GENmu_pt);
                onia_tree->Branch("GENmu_eta",&GENmu_eta);
                onia_tree->Branch("GENmu_phi",&GENmu_phi);
                onia_tree->Branch("GENmu_mass",&GENmu_mass);
                onia_tree->Branch("GENmu_id",&GENmu_id);
                onia_tree->Branch("GENmu_status",&GENmu_status);
                onia_tree->Branch("GENmu_MomId",&GENmu_MomId);
                onia_tree->Branch("GENmu_MomMomId",&GENmu_MomMomId);
                onia_tree->Branch("GENups_DaughtersId",&GENups_DaughtersId);
                onia_tree->Branch("GENups_Daughter_mupt",&GENups_Daughter_mupt);
                onia_tree->Branch("GENups_Daughter_mueta",&GENups_Daughter_mueta);
                onia_tree->Branch("GENups_Daughter_muphi",&GENups_Daughter_muphi);
                onia_tree->Branch("GENups_Daughter_mumass",&GENups_Daughter_mumass);
                onia_tree->Branch("GENups_Daughter_mustatus",&GENups_Daughter_mustatus);
                onia_tree->Branch("GENups_MomId",&GENups_MomId);
                onia_tree->Branch("GENups_pt",&GENups_pt);
                onia_tree->Branch("GENdimu_mass",&GENdimu_mass);
                onia_tree->Branch("GENdimu_pt",&GENdimu_pt);
                onia_tree->Branch("GENdimu_eta",&GENdimu_eta);
                onia_tree->Branch("GENdimu_y",&GENdimu_y);
                onia_tree->Branch("GENdimu_phi",&GENdimu_phi);
                onia_tree->Branch("GENups_eta",&GENups_eta);
                onia_tree->Branch("GENups_y",&GENups_y);
                onia_tree->Branch("GENups_phi",&GENups_phi);
                onia_tree->Branch("GENups_mass",&GENups_mass);
                onia_tree->Branch("GENX_mass",&GENX_mass);
                onia_tree->Branch("GENX_pt",&GENX_pt);
                onia_tree->Branch("GENX_eta",&GENX_eta);
                onia_tree->Branch("GENX_y",&GENX_y);
                onia_tree->Branch("GENX_phi",&GENX_phi);
		onia_tree->Branch("fourMuFit_Mass_allComb",&fourMuFit_Mass_allComb);
		onia_tree->Branch("fourMuFit_Mass",&fourMuFit_Mass);
		onia_tree->Branch("fourMuFit_MassErr",&fourMuFit_MassErr);
		onia_tree->Branch("fourMuFit_Pt",&fourMuFit_Pt);
		onia_tree->Branch("fourMuFit_Eta",&fourMuFit_Eta);
		onia_tree->Branch("fourMuFit_Phi",&fourMuFit_Phi);
		onia_tree->Branch("fourMuFit_VtxX",&fourMuFit_VtxX);
		onia_tree->Branch("fourMuFit_VtxY",&fourMuFit_VtxY);
		onia_tree->Branch("fourMuFit_VtxZ",&fourMuFit_VtxZ);
                onia_tree->Branch("fourMuFit_VtxXE",&fourMuFit_VtxXE);
                onia_tree->Branch("fourMuFit_VtxYE",&fourMuFit_VtxYE);
                onia_tree->Branch("fourMuFit_VtxZE",&fourMuFit_VtxZE);
                onia_tree->Branch("fourMuFit_PVX",&fourMuFit_PVX);
                onia_tree->Branch("fourMuFit_PVY",&fourMuFit_PVY);
                onia_tree->Branch("fourMuFit_PVZ",&fourMuFit_PVZ);
                onia_tree->Branch("fourMuFit_PVXE",&fourMuFit_PVXE);
                onia_tree->Branch("fourMuFit_PVYE",&fourMuFit_PVYE);
                onia_tree->Branch("fourMuFit_PVZE",&fourMuFit_PVZE);
                onia_tree->Branch("fourMuFit_PV_XYE",&fourMuFit_PV_XYE);
                onia_tree->Branch("fourMuFit_PV_XZE",&fourMuFit_PV_XZE);
                onia_tree->Branch("fourMuFit_PV_YZE",&fourMuFit_PV_YZE);
                onia_tree->Branch("fourMuFit_PV_CL",&fourMuFit_PV_CL);
                onia_tree->Branch("fourMuFit_VtxProb",&fourMuFit_VtxProb);

                //Jinfeng_Xdistance
                onia_tree->Branch("fourMuFit_cosAlpha_noMC", &fourMuFit_cosAlpha_noMC);
                onia_tree->Branch("fourMuFit_lxy_noMC", &fourMuFit_lxy_noMC);
                onia_tree->Branch("fourMuFit_lxyerr_noMC", &fourMuFit_lxyerr_noMC);
                onia_tree->Branch("fourMuFit_lxySig_noMC", &fourMuFit_lxySig_noMC);
                onia_tree->Branch("fourMuFit_ctauPV_noMC", &fourMuFit_ctauPV_noMC);
                onia_tree->Branch("fourMuFit_ctauErrPV_noMC", &fourMuFit_ctauErrPV_noMC);
                onia_tree->Branch("fourMuFit_lxyPV_noMC", &fourMuFit_lxyPV_noMC);
                onia_tree->Branch("fourMuFit_cosAlpha", &fourMuFit_cosAlpha);
                onia_tree->Branch("fourMuFit_lxy", &fourMuFit_lxy);
                onia_tree->Branch("fourMuFit_lxyerr", &fourMuFit_lxyerr);
                onia_tree->Branch("fourMuFit_lxySig", &fourMuFit_lxySig);
                onia_tree->Branch("fourMuFit_ctauPV", &fourMuFit_ctauPV);
                onia_tree->Branch("fourMuFit_ctauErrPV", &fourMuFit_ctauErrPV);
                onia_tree->Branch("fourMuFit_lxyPV", &fourMuFit_lxyPV);

                onia_tree->Branch("fourMuFit_Double_ups_VtxProb",&fourMuFit_Double_ups_VtxProb);
                onia_tree->Branch("fourMuFit_ups1_VtxProb",&fourMuFit_ups1_VtxProb);
                onia_tree->Branch("fourMuFit_ups2_VtxProb",&fourMuFit_ups2_VtxProb);
                onia_tree->Branch("fourMuFit_ups1_VtxProb_noMC",&fourMuFit_ups1_VtxProb_noMC);
                onia_tree->Branch("fourMuFit_ups2_VtxProb_noMC",&fourMuFit_ups2_VtxProb_noMC);
                onia_tree->Branch("fourMuFit_ups1_pt",&fourMuFit_ups1_pt);
                onia_tree->Branch("fourMuFit_ups2_pt",&fourMuFit_ups2_pt);
                onia_tree->Branch("fourMuFit_ups1_rapidity",&fourMuFit_ups1_rapidity);
                onia_tree->Branch("fourMuFit_ups2_rapidity",&fourMuFit_ups2_rapidity);
                onia_tree->Branch("fourMuFit_ups1_mass",&fourMuFit_ups1_mass);
                onia_tree->Branch("fourMuFit_ups2_mass",&fourMuFit_ups2_mass);
                onia_tree->Branch("fourMuFit_ups1_massError",&fourMuFit_ups1_massError);
                onia_tree->Branch("fourMuFit_ups2_massError",&fourMuFit_ups2_massError);
                onia_tree->Branch("fourMuFit_ups1_cTau_noMC",&fourMuFit_ups1_cTau_noMC);
                onia_tree->Branch("fourMuFit_ups2_cTau_noMC",&fourMuFit_ups2_cTau_noMC);
                onia_tree->Branch("fourMuFit_ups1_cTau_MC",&fourMuFit_ups1_cTau_MC);
                onia_tree->Branch("fourMuFit_ups2_cTau_MC",&fourMuFit_ups2_cTau_MC);
                onia_tree->Branch("fourMuFit_ups1_cTauErr_noMC",&fourMuFit_ups1_cTauErr_noMC);
                onia_tree->Branch("fourMuFit_ups2_cTauErr_noMC",&fourMuFit_ups2_cTauErr_noMC);
                onia_tree->Branch("fourMuFit_ups1_cTauErr_MC",&fourMuFit_ups1_cTauErr_MC);
                onia_tree->Branch("fourMuFit_ups2_cTauErr_MC",&fourMuFit_ups2_cTauErr_MC);

                //Jinfeng
                onia_tree->Branch("fourMuFit_ups1_Lxy_noMC", &fourMuFit_ups1_Lxy_noMC);
                onia_tree->Branch("fourMuFit_ups2_Lxy_noMC", &fourMuFit_ups2_Lxy_noMC);
                onia_tree->Branch("fourMuFit_ups1_Lxy_MC", &fourMuFit_ups1_Lxy_MC);
                onia_tree->Branch("fourMuFit_ups2_Lxy_MC", &fourMuFit_ups2_Lxy_MC);
               
                onia_tree->Branch("fourMuFit_ups1_CosAlpha_noMC", &fourMuFit_ups1_CosAlpha_noMC);
                onia_tree->Branch("fourMuFit_ups2_CosAlpha_noMC", &fourMuFit_ups2_CosAlpha_noMC);
                onia_tree->Branch("fourMuFit_ups1_CosAlpha_MC", &fourMuFit_ups1_CosAlpha_MC);
                onia_tree->Branch("fourMuFit_ups2_CosAlpha_MC", &fourMuFit_ups2_CosAlpha_MC);

                onia_tree->Branch("fourMuFit_ups1_LxyPV_noMC", &fourMuFit_ups1_LxyPV_noMC);
                onia_tree->Branch("fourMuFit_ups2_LxyPV_noMC", &fourMuFit_ups2_LxyPV_noMC);
                onia_tree->Branch("fourMuFit_ups1_LxyPV_MC", &fourMuFit_ups1_LxyPV_MC);
                onia_tree->Branch("fourMuFit_ups2_LxyPV_MC", &fourMuFit_ups2_LxyPV_MC);

                //Jinfeng 10.3
                onia_tree->Branch("fourMuFit_ups1_LxyPVErr_noMC", &fourMuFit_ups1_LxyPVErr_noMC);
                onia_tree->Branch("fourMuFit_ups2_LxyPVErr_noMC", &fourMuFit_ups2_LxyPVErr_noMC);
                onia_tree->Branch("fourMuFit_ups1_LxyPVErr_MC", &fourMuFit_ups1_LxyPVErr_MC);
                onia_tree->Branch("fourMuFit_ups2_LxyPVErr_MC", &fourMuFit_ups2_LxyPVErr_MC);
                onia_tree->Branch("fourMuFit_ups1_LxyPVSig_noMC", &fourMuFit_ups1_LxyPVSig_noMC);
                onia_tree->Branch("fourMuFit_ups2_LxyPVSig_noMC", &fourMuFit_ups2_LxyPVSig_noMC);
                onia_tree->Branch("fourMuFit_ups1_LxyPVSig_MC", &fourMuFit_ups1_LxyPVSig_MC);
                onia_tree->Branch("fourMuFit_ups2_LxyPVSig_MC", &fourMuFit_ups2_LxyPVSig_MC);
                onia_tree->Branch("fourMuFit_Distance_noMC", &fourMuFit_Distance_noMC);
                onia_tree->Branch("fourMuFit_Distance_MC", &fourMuFit_Distance_MC);
                onia_tree->Branch("fourMuFit_DistanceErr_noMC", &fourMuFit_DistanceErr_noMC);
                onia_tree->Branch("fourMuFit_DistanceErr_MC", &fourMuFit_DistanceErr_MC);
                onia_tree->Branch("fourMuFit_DistanceSig_noMC", &fourMuFit_DistanceSig_noMC);
                onia_tree->Branch("fourMuFit_DistanceSig_MC", &fourMuFit_DistanceSig_MC);

                //Jinfeng 12.14
                onia_tree->Branch("Trigger_Matched", &Trigger_Matched);

                onia_tree->Branch("fourMuFit_ups1_VtxX",&fourMuFit_ups1_VtxX);
                onia_tree->Branch("fourMuFit_ups1_VtxY",&fourMuFit_ups1_VtxY);
                onia_tree->Branch("fourMuFit_ups1_VtxZ",&fourMuFit_ups1_VtxZ);
                onia_tree->Branch("fourMuFit_ups1_VtxXE",&fourMuFit_ups1_VtxXE);
                onia_tree->Branch("fourMuFit_ups1_VtxYE",&fourMuFit_ups1_VtxYE);
                onia_tree->Branch("fourMuFit_ups1_VtxZE",&fourMuFit_ups1_VtxZE);
                onia_tree->Branch("fourMuFit_ups2_VtxX",&fourMuFit_ups2_VtxX);
                onia_tree->Branch("fourMuFit_ups2_VtxY",&fourMuFit_ups2_VtxY);
                onia_tree->Branch("fourMuFit_ups2_VtxZ",&fourMuFit_ups2_VtxZ);
                onia_tree->Branch("fourMuFit_ups2_VtxXE",&fourMuFit_ups2_VtxXE);
                onia_tree->Branch("fourMuFit_ups2_VtxYE",&fourMuFit_ups2_VtxYE);
                onia_tree->Branch("fourMuFit_ups2_VtxZE",&fourMuFit_ups2_VtxZE);
                onia_tree->Branch("fourMuFit_wrong_ups1_mass",&fourMuFit_wrong_ups1_mass);
                onia_tree->Branch("fourMuFit_wrong_ups2_mass",&fourMuFit_wrong_ups2_mass);
                onia_tree->Branch("fourMuFit_wrong_ups1_massError",&fourMuFit_wrong_ups1_massError);
                onia_tree->Branch("fourMuFit_wrong_ups2_massError",&fourMuFit_wrong_ups2_massError);
		onia_tree->Branch("fourMuFit_Chi2",&fourMuFit_Chi2);
		onia_tree->Branch("fourMuFit_ndof",&fourMuFit_ndof);
                onia_tree->Branch("fourMuFit_Double_ups_Chi2",&fourMuFit_Double_ups_Chi2);
                onia_tree->Branch("fourMuFit_Double_ups_ndof",&fourMuFit_Double_ups_ndof);            
		onia_tree->Branch("fourMuFit_mu1Pt",&fourMuFit_mu1Pt);
		onia_tree->Branch("fourMuFit_mu1Eta",&fourMuFit_mu1Eta);
		onia_tree->Branch("fourMuFit_mu1Phi",&fourMuFit_mu1Phi);
		onia_tree->Branch("fourMuFit_mu1E",&fourMuFit_mu1E);
                onia_tree->Branch("fourMuFit_mu1frompv",&fourMuFit_mu1frompv);
		onia_tree->Branch("fourMuFit_mu2Pt",&fourMuFit_mu2Pt);
		onia_tree->Branch("fourMuFit_mu2Eta",&fourMuFit_mu2Eta);
		onia_tree->Branch("fourMuFit_mu2Phi",&fourMuFit_mu2Phi);
		onia_tree->Branch("fourMuFit_mu2E",&fourMuFit_mu2E);
                onia_tree->Branch("fourMuFit_mu2frompv",&fourMuFit_mu2frompv);
		onia_tree->Branch("fourMuFit_mu3Pt",&fourMuFit_mu3Pt);
		onia_tree->Branch("fourMuFit_mu3Eta",&fourMuFit_mu3Eta);
		onia_tree->Branch("fourMuFit_mu3Phi",&fourMuFit_mu3Phi);
		onia_tree->Branch("fourMuFit_mu3E",&fourMuFit_mu3E);
                onia_tree->Branch("fourMuFit_mu3frompv",&fourMuFit_mu3frompv);
		onia_tree->Branch("fourMuFit_mu4Pt",&fourMuFit_mu4Pt);
		onia_tree->Branch("fourMuFit_mu4Eta",&fourMuFit_mu4Eta);
		onia_tree->Branch("fourMuFit_mu4Phi",&fourMuFit_mu4Phi);
		onia_tree->Branch("fourMuFit_mu4E",&fourMuFit_mu4E);
                onia_tree->Branch("fourMuFit_mu4frompv",&fourMuFit_mu4frompv);
                onia_tree->Branch("fourMuFit_mu3_trg_dR",&fourMuFit_mu3_trg_dR);
                onia_tree->Branch("fourMuFit_mu4_trg_dR",&fourMuFit_mu4_trg_dR);
                onia_tree->Branch("AllRecoMuons_Pt",   &AllRecoMuons_Pt);
                onia_tree->Branch("AllRecoMuons_Eta",   &AllRecoMuons_Eta);
                onia_tree->Branch("mu1_Pt",   &mu1_Pt);
                onia_tree->Branch("mu1_Eta",  &mu1_Eta);
                onia_tree->Branch("mu1_Phi",  &mu1_Phi);
                onia_tree->Branch("mu1_E",    &mu1_E);
                onia_tree->Branch("mu1_genindex",   &mu1_genindex);
                onia_tree->Branch("mu2_Pt",   &mu2_Pt);
                onia_tree->Branch("mu2_Eta",  &mu2_Eta);
                onia_tree->Branch("mu2_Phi",  &mu2_Phi);
                onia_tree->Branch("mu2_E",    &mu2_E);
                onia_tree->Branch("mu2_genindex",   &mu2_genindex);
		onia_tree->Branch("mu3_Pt",   &mu3_Pt);
		onia_tree->Branch("mu3_Eta",  &mu3_Eta);
		onia_tree->Branch("mu3_Phi",  &mu3_Phi);
		onia_tree->Branch("mu3_E",    &mu3_E);
                //onia_tree->Branch("mu3_genindex",   &mu3_genindex);
		onia_tree->Branch("mu4_Pt",   &mu4_Pt);
		onia_tree->Branch("mu4_Eta",  &mu4_Eta);
		onia_tree->Branch("mu4_Phi",  &mu4_Phi);
		onia_tree->Branch("mu4_E",    &mu4_E);
                //onia_tree->Branch("mu4_genindex",   &mu4_genindex);
                onia_tree->Branch("mu1Charge", &mu1Charge);
                onia_tree->Branch("mu2Charge",   &mu2Charge);
		onia_tree->Branch("mu3Charge", &mu3Charge);
		onia_tree->Branch("mu4Charge",   &mu4Charge);
		onia_tree->Branch("mu3_d0",   &mu3_d0);
		onia_tree->Branch("mu3_d0err",   &mu3_d0err);
		onia_tree->Branch("mu4_d0",   &mu4_d0);
		onia_tree->Branch("mu4_d0err",   &mu4_d0err);
		onia_tree->Branch("mu3_dz",   &mu3_dz);
		onia_tree->Branch("mu3_dzerr",   &mu3_dzerr);
		onia_tree->Branch("mu4_dz",   &mu4_dz);
		onia_tree->Branch("mu4_dzerr",   &mu4_dzerr);
                onia_tree->Branch("mu1_Soft",   &mu1_Soft);
                onia_tree->Branch("mu2_Soft",   &mu2_Soft);
                onia_tree->Branch("mu3_Soft",   &mu3_Soft);
                onia_tree->Branch("mu4_Soft",   &mu4_Soft);
		onia_tree->Branch("mu1_Tight",   &mu1_Tight);
		onia_tree->Branch("mu2_Tight",   &mu2_Tight);
		onia_tree->Branch("mu3_Tight",   &mu3_Tight);
		onia_tree->Branch("mu4_Tight",   &mu4_Tight);
		onia_tree->Branch("mu1_Medium",   &mu1_Medium);
		onia_tree->Branch("mu2_Medium",   &mu2_Medium);
		onia_tree->Branch("mu3_Medium",   &mu3_Medium);
		onia_tree->Branch("mu4_Medium",   &mu4_Medium);
		onia_tree->Branch("mu1_Loose",   &mu1_Loose);
		onia_tree->Branch("mu2_Loose",   &mu2_Loose);
		onia_tree->Branch("mu3_Loose",   &mu3_Loose);
		onia_tree->Branch("mu4_Loose",   &mu4_Loose);
		onia_tree->Branch("mu1_pdgID",   &mu1_pdgID);
		onia_tree->Branch("mu2_pdgID",   &mu2_pdgID);
		onia_tree->Branch("mu3_pdgID",   &mu3_pdgID);
		onia_tree->Branch("mu4_pdgID",   &mu4_pdgID);
                onia_tree->Branch("fourMuFit_mu12overlap", &fourMuFit_mu12overlap);
                onia_tree->Branch("fourMuFit_mu13overlap", &fourMuFit_mu13overlap);
                onia_tree->Branch("fourMuFit_mu14overlap", &fourMuFit_mu14overlap);
                onia_tree->Branch("fourMuFit_mu23overlap", &fourMuFit_mu23overlap);
                onia_tree->Branch("fourMuFit_mu24overlap", &fourMuFit_mu24overlap);
                onia_tree->Branch("fourMuFit_mu34overlap", &fourMuFit_mu34overlap);
                onia_tree->Branch("fourMuFit_mu12SharedSeg", &fourMuFit_mu12SharedSeg);
                onia_tree->Branch("fourMuFit_mu13SharedSeg", &fourMuFit_mu13SharedSeg);
                onia_tree->Branch("fourMuFit_mu14SharedSeg", &fourMuFit_mu14SharedSeg);
                onia_tree->Branch("fourMuFit_mu23SharedSeg", &fourMuFit_mu23SharedSeg);
                onia_tree->Branch("fourMuFit_mu24SharedSeg", &fourMuFit_mu24SharedSeg);
                onia_tree->Branch("fourMuFit_mu34SharedSeg", &fourMuFit_mu34SharedSeg);
		/*
			onia_tree->Branch("fourMuFit_Mass_allComb",&fourMuFit_Mass_allComb);
			onia_tree->Branch("fourMuFit_Mass",&fourMuFit_Mass,"fourMuFit_Mass/F");
			onia_tree->Branch("fourMuFit_MassErr",&fourMuFit_MassErr,"fourMuFit_MassErr/F");
			onia_tree->Branch("fourMuFit_VtxX",&fourMuFit_VtxX,"fourMuFit_VtxX/F");
			onia_tree->Branch("fourMuFit_VtxY",&fourMuFit_VtxY,"fourMuFit_VtxY/F");
			onia_tree->Branch("fourMuFit_VtxZ",&fourMuFit_VtxZ,"fourMuFit_VtxZ/F");
			onia_tree->Branch("fourMuFit_VtxProb",&fourMuFit_VtxProb,"fourMuFit_VtxProb/F");
			onia_tree->Branch("fourMuFit_Chi2",&fourMuFit_Chi2,"fourMuFit_Chi2/F");
			onia_tree->Branch("fourMuFit_ndof",&fourMuFit_ndof,"fourMuFit_ndof/I");
			onia_tree->Branch("fourMuFit_p4",  "TLorentzVector", &fourMuFit_p4);
			onia_tree->Branch("fourMuFit_mu1p4",  "TLorentzVector", &fourMuFit_mu1p4);
			onia_tree->Branch("fourMuFit_mu2p4",  "TLorentzVector", &fourMuFit_mu2p4);
			onia_tree->Branch("fourMuFit_mu3p4",  "TLorentzVector", &fourMuFit_mu3p4);
			onia_tree->Branch("fourMuFit_mu4p4",  "TLorentzVector", &fourMuFit_mu4p4);
			onia_tree->Branch("mu3Charge",   &mu3Charge,    "mu3Charge/I");
			onia_tree->Branch("mu4Charge",   &mu4Charge,    "mu4Charge/I");
			onia_tree->Branch("mu3_p4",  "TLorentzVector", &mu3_p4);
			onia_tree->Branch("mu4_p4",  "TLorentzVector", &mu4_p4);
			onia_tree->Branch("mu3_d0",   &mu3_d0,    "mu3_d0/F");
			onia_tree->Branch("mu3_d0err",   &mu3_d0err,    "mu3_d0err/F");
			onia_tree->Branch("mu4_d0",   &mu4_d0,    "mu4_d0/F");
			onia_tree->Branch("mu4_d0err",   &mu4_d0err,    "mu4_d0err/F");
			onia_tree->Branch("mu3_dz",   &mu3_dz,    "mu3_dz/F");
			onia_tree->Branch("mu3_dzerr",   &mu3_dzerr,    "mu3_dzerr/F");
			onia_tree->Branch("mu4_dz",   &mu4_dz,    "mu4_dz/F");
			onia_tree->Branch("mu4_dzerr",   &mu4_dzerr,    "mu4_dzerr/F");
			onia_tree->Branch("mu1_Tight",   &mu1_Tight,    "mu1_Tight/I");
			onia_tree->Branch("mu2_Tight",   &mu2_Tight,    "mu2_Tight/I");
			onia_tree->Branch("mu3_Tight",   &mu3_Tight,    "mu3_Tight/I");
			onia_tree->Branch("mu4_Tight",   &mu4_Tight,    "mu4_Tight/I");
			onia_tree->Branch("mu3_pdgID",   &mu3_pdgID,    "mu3_pdgID/I");
			onia_tree->Branch("mu4_pdgID",   &mu4_pdgID,    "mu4_pdgID/I");
			onia_tree->Branch("upsilonMuons", &upsilonMuons);	
			onia_tree->Branch("theRestMuons", &theRestMuons); 
			*/
		onia_tree->Branch("mu1_p4_bestYMass",  "TLorentzVector", &mu1_p4_bestYMass);
		onia_tree->Branch("mu2_p4_bestYMass",  "TLorentzVector", &mu2_p4_bestYMass);
		onia_tree->Branch("mu1Charge_bestYMass",   &mu1Charge_bestYMass,    "mu1Charge_bestYMass/I");
		onia_tree->Branch("mu2Charge_bestYMass",   &mu2Charge_bestYMass,    "mu2Charge_bestYMass/I");
		onia_tree->Branch("mu1_d0_bestYMass",   &mu1_d0_bestYMass,    "mu1_d0_bestYMass/F");
		onia_tree->Branch("mu1_d0err_bestYMass",   &mu1_d0err_bestYMass,    "mu1_d0err_bestYMass/F");
		onia_tree->Branch("mu2_d0_bestYMass",   &mu2_d0_bestYMass,    "mu2_d0_bestYMass/F");
		onia_tree->Branch("mu2_d0err_bestYMass",   &mu2_d0err_bestYMass,    "mu2_d0err_bestYMass/F");
		onia_tree->Branch("mu1_dz_bestYMass",   &mu1_dz_bestYMass,    "mu1_dz_bestYMass/F");
		onia_tree->Branch("mu1_dzerr_bestYMass",   &mu1_dzerr_bestYMass,    "mu1_dzerr_bestYMass/F");
		onia_tree->Branch("mu2_dz_bestYMass",   &mu2_dz_bestYMass,    "mu2_dz_bestYMass/F");
		onia_tree->Branch("mu2_dzerr_bestYMass",   &mu2_dzerr_bestYMass,    "mu2_dzerr_bestYMass/F");
		onia_tree->Branch("dimuon_p4_bestYMass", "TLorentzVector", &dimuon_p4_bestYMass);
		onia_tree->Branch("mumufit_Mass_bestYMass",&mumufit_Mass_bestYMass,"mumufit_Mass_bestYMass/F");
		onia_tree->Branch("mumufit_MassErr_bestYMass",&mumufit_MassErr_bestYMass,"mumufit_MassErr_bestYMass/F");
		onia_tree->Branch("mumufit_VtxCL_bestYMass",&mumufit_VtxCL_bestYMass,"mumufit_VtxCL_bestYMass/F");
		onia_tree->Branch("mumufit_VtxCL2_bestYMass",&mumufit_VtxCL2_bestYMass,"mumufit_VtxCL2_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxX_bestYMass",&mumufit_DecayVtxX_bestYMass,"mumufit_DecayVtxX_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxY_bestYMass",&mumufit_DecayVtxY_bestYMass,"mumufit_DecayVtxY_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxZ_bestYMass",&mumufit_DecayVtxZ_bestYMass,"mumufit_DecayVtxZ_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxXE_bestYMass",&mumufit_DecayVtxXE_bestYMass,"mumufit_DecayVtxXE_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxYE_bestYMass",&mumufit_DecayVtxYE_bestYMass,"mumufit_DecayVtxYE_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxZE_bestYMass",&mumufit_DecayVtxZE_bestYMass,"mumufit_DecayVtxZE_bestYMass/F");
		onia_tree->Branch("mumufit_p4_bestYMass",  "TLorentzVector", &mumufit_p4_bestYMass);
		onia_tree->Branch("bestVertex_and_bestYMass", &bestVertex_and_bestYMass,"bestVertex_and_bestYMass/I");

		onia_tree->Branch("fourMuFit_Mass_allComb_mix_bestYMass",&fourMuFit_Mass_allComb_mix_bestYMass);
		onia_tree->Branch("fourMuFit_Mass_mix_bestYMass",&fourMuFit_Mass_mix_bestYMass,"fourMuFit_Mass_mix_bestYMass/F");
		onia_tree->Branch("fourMuFit_MassErr_mix_bestYMass",&fourMuFit_MassErr_mix_bestYMass,"fourMuFit_MassErr_mix_bestYMass/F");
		onia_tree->Branch("fourMuFit_VtxX_mix_bestYMass",&fourMuFit_VtxX_mix_bestYMass,"fourMuFit_VtxX_mix_bestYMass/F");
		onia_tree->Branch("fourMuFit_VtxY_mix_bestYMass",&fourMuFit_VtxY_mix_bestYMass,"fourMuFit_VtxY_mix_bestYMass/F");
		onia_tree->Branch("fourMuFit_VtxZ_mix_bestYMass",&fourMuFit_VtxZ_mix_bestYMass,"fourMuFit_VtxZ_mix_bestYMass/F");
		onia_tree->Branch("fourMuFit_VtxProb_mix_bestYMass",&fourMuFit_VtxProb_mix_bestYMass,"fourMuFit_VtxProb_mix_bestYMass/F");
		onia_tree->Branch("fourMuFit_Chi2_mix_bestYMass",&fourMuFit_Chi2_mix_bestYMass,"fourMuFit_Chi2_mix_bestYMass/F");
		onia_tree->Branch("fourMuFit_ndof_mix_bestYMass",&fourMuFit_ndof_mix_bestYMass,"fourMuFit_ndof_mix_bestYMass/I");
		onia_tree->Branch("fourMuFit_3plus1_mix_bestYMass",&fourMuFit_3plus1_mix_bestYMass,"fourMuFit_3plus1_mix_bestYMass/I");
		onia_tree->Branch("fourMuFit_p4_mix_bestYMass",  "TLorentzVector", &fourMuFit_p4_mix_bestYMass);
		onia_tree->Branch("fourMuFit_mu1p4_mix_bestYMass",  "TLorentzVector", &fourMuFit_mu1p4_mix_bestYMass);
		onia_tree->Branch("fourMuFit_mu2p4_mix_bestYMass",  "TLorentzVector", &fourMuFit_mu2p4_mix_bestYMass);
		onia_tree->Branch("fourMuFit_mu3p4_mix_bestYMass",  "TLorentzVector", &fourMuFit_mu3p4_mix_bestYMass);
		onia_tree->Branch("fourMuFit_mu4p4_mix_bestYMass",  "TLorentzVector", &fourMuFit_mu4p4_mix_bestYMass);
		onia_tree->Branch("mu3Charge_mix_bestYMass",   &mu3Charge_mix_bestYMass,    "mu3Charge_mix_bestYMass/I");
		onia_tree->Branch("mu4Charge_mix_bestYMass",   &mu4Charge_mix_bestYMass,    "mu4Charge_mix_bestYMass/I");
		onia_tree->Branch("mu3_p4_mix_bestYMass",  "TLorentzVector", &mu3_p4_mix_bestYMass);
		onia_tree->Branch("mu4_p4_mix_bestYMass",  "TLorentzVector", &mu4_p4_mix_bestYMass);
		onia_tree->Branch("mu3_d0_mix_bestYMass",   &mu3_d0_mix_bestYMass,    "mu3_d0_mix_bestYMass/F");
		onia_tree->Branch("mu3_d0err_mix_bestYMass",   &mu3_d0err_mix_bestYMass,    "mu3_d0err_mix_bestYMass/F");
		onia_tree->Branch("mu4_d0_mix_bestYMass",   &mu4_d0_mix_bestYMass,    "mu4_d0_mix_bestYMass/F");
		onia_tree->Branch("mu4_d0err_mix_bestYMass",   &mu4_d0err_mix_bestYMass,    "mu4_d0err_mix_bestYMass/F");
		onia_tree->Branch("mu3_dz_mix_bestYMass",   &mu3_dz_mix_bestYMass,    "mu3_dz_mix_bestYMass/F");
		onia_tree->Branch("mu3_dzerr_mix_bestYMass",   &mu3_dzerr_mix_bestYMass,    "mu3_dzerr_mix_bestYMass/F");
		onia_tree->Branch("mu4_dz_mix_bestYMass",   &mu4_dz_mix_bestYMass,    "mu4_dz_mix_bestYMass/F");
		onia_tree->Branch("mu4_dzerr_mix_bestYMass",   &mu4_dzerr_mix_bestYMass,    "mu4_dzerr_mix_bestYMass/F");

		onia_tree->Branch("fourMuFit_Mass_allComb_bestYMass",&fourMuFit_Mass_allComb_bestYMass);
		onia_tree->Branch("fourMuFit_Mass_bestYMass",&fourMuFit_Mass_bestYMass,"fourMuFit_Mass_bestYMass/F");
		onia_tree->Branch("fourMuFit_MassErr_bestYMass",&fourMuFit_MassErr_bestYMass,"fourMuFit_MassErr_bestYMass/F");
		onia_tree->Branch("fourMuFit_VtxX_bestYMass",&fourMuFit_VtxX_bestYMass,"fourMuFit_VtxX_bestYMass/F");
		onia_tree->Branch("fourMuFit_VtxY_bestYMass",&fourMuFit_VtxY_bestYMass,"fourMuFit_VtxY_bestYMass/F");
		onia_tree->Branch("fourMuFit_VtxZ_bestYMass",&fourMuFit_VtxZ_bestYMass,"fourMuFit_VtxZ_bestYMass/F");
		onia_tree->Branch("fourMuFit_VtxProb_bestYMass",&fourMuFit_VtxProb_bestYMass,"fourMuFit_VtxProb_bestYMass/F");
		onia_tree->Branch("fourMuFit_Chi2_bestYMass",&fourMuFit_Chi2_bestYMass,"fourMuFit_Chi2_bestYMass/F");
		onia_tree->Branch("fourMuFit_ndof_bestYMass",&fourMuFit_ndof_bestYMass,"fourMuFit_ndof_bestYMass/I");
		onia_tree->Branch("fourMuFit_p4_bestYMass",  "TLorentzVector", &fourMuFit_p4_bestYMass);
		onia_tree->Branch("fourMuFit_mu1p4_bestYMass",  "TLorentzVector", &fourMuFit_mu1p4_bestYMass);
		onia_tree->Branch("fourMuFit_mu2p4_bestYMass",  "TLorentzVector", &fourMuFit_mu2p4_bestYMass);
		onia_tree->Branch("fourMuFit_mu3p4_bestYMass",  "TLorentzVector", &fourMuFit_mu3p4_bestYMass);
		onia_tree->Branch("fourMuFit_mu4p4_bestYMass",  "TLorentzVector", &fourMuFit_mu4p4_bestYMass);
		onia_tree->Branch("mu3Charge_bestYMass",   &mu3Charge_bestYMass,    "mu3Charge_bestYMass/I");
		onia_tree->Branch("mu4Charge_bestYMass",   &mu4Charge_bestYMass,    "mu4Charge_bestYMass/I");
		onia_tree->Branch("mu3_p4_bestYMass",  "TLorentzVector", &mu3_p4_bestYMass);
		onia_tree->Branch("mu4_p4_bestYMass",  "TLorentzVector", &mu4_p4_bestYMass);
		onia_tree->Branch("mu3_d0_bestYMass",   &mu3_d0_bestYMass,    "mu3_d0_bestYMass/F");
		onia_tree->Branch("mu3_d0err_bestYMass",   &mu3_d0err_bestYMass,    "mu3_d0err_bestYMass/F");
		onia_tree->Branch("mu4_d0_bestYMass",   &mu4_d0_bestYMass,    "mu4_d0_bestYMass/F");
		onia_tree->Branch("mu4_d0err_bestYMass",   &mu4_d0err_bestYMass,    "mu4_d0err_bestYMass/F");
		onia_tree->Branch("mu3_dz_bestYMass",   &mu3_dz_bestYMass,    "mu3_dz_bestYMass/F");
		onia_tree->Branch("mu3_dzerr_bestYMass",   &mu3_dzerr_bestYMass,    "mu3_dzerr_bestYMass/F");
		onia_tree->Branch("mu4_dz_bestYMass",   &mu4_dz_bestYMass,    "mu4_dz_bestYMass/F");
		onia_tree->Branch("mu4_dzerr_bestYMass",   &mu4_dzerr_bestYMass,    "mu4_dzerr_bestYMass/F");
		onia_tree->Branch("mu1_Tight_bestYMass",   &mu1_Tight_bestYMass,    "mu1_Tight_bestYMass/I");
		onia_tree->Branch("mu2_Tight_bestYMass",   &mu2_Tight_bestYMass,    "mu2_Tight_bestYMass/I");
		onia_tree->Branch("mu3_Tight_bestYMass",   &mu3_Tight_bestYMass,    "mu3_Tight_bestYMass/I");
		onia_tree->Branch("mu4_Tight_bestYMass",   &mu4_Tight_bestYMass,    "mu4_Tight_bestYMass/I");
		onia_tree->Branch("mu3_pdgID_bestYMass",   &mu3_pdgID_bestYMass,    "mu3_pdgID_bestYMass/I");
		onia_tree->Branch("mu4_pdgID_bestYMass",   &mu4_pdgID_bestYMass,    "mu4_pdgID_bestYMass/I");
                onia_tree->Branch("mu1_Medium_bestYMass",&mu1_Medium_bestYMass,"mu1_Medium_bestYMass/I");
                onia_tree->Branch("mu2_Medium_bestYMass",&mu2_Medium_bestYMass,"mu2_Medium_bestYMass/I");
                onia_tree->Branch("mu3_Medium_bestYMass",&mu3_Medium_bestYMass,"mu3_Medium_bestYMass/I");
                onia_tree->Branch("mu4_Medium_bestYMass",&mu4_Medium_bestYMass,"mu4_Medium_bestYMass/I");
                onia_tree->Branch("mu1_Loose_bestYMass",&mu1_Loose_bestYMass,"mu1_Loose_bestYMass/I");
                onia_tree->Branch("mu2_Loose_bestYMass",&mu2_Loose_bestYMass,"mu2_Loose_bestYMass/I");
                onia_tree->Branch("mu3_Loose_bestYMass",&mu3_Loose_bestYMass,"mu3_Loose_bestYMass/I");
                onia_tree->Branch("mu4_Loose_bestYMass",&mu4_Loose_bestYMass,"mu4_Loose_bestYMass/I");

	}

	if (isMC_ || OnlyGen_) {
		onia_tree->Branch("mother_pdgId",  &mother_pdgId,     "mother_pdgId/I");
		onia_tree->Branch("dimuon_pdgId",  &dimuon_pdgId,     "dimuon_pdgId/I");
		onia_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
		onia_tree->Branch("gen_mu1_p4",  "TLorentzVector",  &gen_mu1_p4);
		onia_tree->Branch("gen_mu2_p4",  "TLorentzVector",  &gen_mu2_p4);
	}
	prunedgenParticlesSrc_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
	//packedgenParticlesSrc_ = consumes<pat::PackedGenParticle>((edm::InputTag)"packedGenParticles");
   //prunedgenParticlesSrc_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("prunedgenParticlesSrc"))),
   // packedgenParticlesSrc_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("packedgenParticlesSrc"))),

        Total_events = fs->make<TH1F>("Total_events","nEvents in Sample",2,0,2);
/*
	myFourMuM_fit = hists.make<TH1F>("myFourMuM_fit","myFourMuM_fit",50000,0,500.0);
	myFourMuM_fit->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}#mu^{+}#mu^{-}} [GeV/c^{2}]");
	myFourMuVtxP_fit = hists.make<TH1F>("myFourMuVtxP_fit","myFourMuVtxP_fit",1000,0,1);
	myFourMuVtxP_fit->GetXaxis()->SetTitle("#mu^{+}#mu^{-}#mu^{+}#mu^{-} vertex fit #chi^{2} probability");
	myDimuonMass_all = hists.make<TH1F>("myDimuonMass_all","myDimuonMass_all",2000,0,20.0);
	myDimuonMass_all->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
	myY1SFitMass_all = hists.make<TH1F>("myY1SFitMass_all","myY1SFitMass_all",2000,0,20.0);
	myY1SFitMass_all->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
*/
}

MuMuGammaRootupler::~MuMuGammaRootupler() {

}

//
// member functions
//

const reco::Candidate* MuMuGammaRootupler::GetAncestor(const reco::Candidate* p) {
	if (p->numberOfMothers()) {
		if  ((p->mother(0))->pdgId() == p->pdgId()) return GetAncestor(p->mother(0));
		else return p->mother(0);
	}
	//std::cout << "GetAncestor: Inconsistet ancestor, particle does not have a mother " << std::endl;
	return p;
}

//Check recursively if any ancestor of particle is the given one
bool MuMuGammaRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
	if (ancestor == particle ) return true;
	for (size_t i=0; i< particle->numberOfMothers(); i++) {
		if (isAncestor(ancestor, particle->mother(i))) return true;
	}
	return false;
}

/* Grab Trigger information. Save it in variable trigger, trigger is an int between 0 and 127, in binary it is:
	(pass 2)(pass 1)(pass 0)
	ex. 7 = pass 0, 1 and 2
	ex. 6 = pass 2, 3
	ex. 1 = pass 0
	*/

bool MuMuGammaRootupler::triggerDecision(edm::Handle<edm::TriggerResults> &hltR, int iTrigger){
  bool triggerPassed = false;
  if(hltR->wasrun(iTrigger) &&
     hltR->accept(iTrigger) &&
     !hltR->error(iTrigger) ){
    triggerPassed = true;
  }
  return triggerPassed;
}


void MuMuGammaRootupler::analyzeTrigger(edm::Handle<edm::TriggerResults> &hltR,
                                       edm::Handle<trigger::TriggerEvent> &hltE,
                                       const std::string& triggerName) {
  using namespace trigger;

  const unsigned int n(hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  if(verbose){
    std::cout<<" n = "<<n<<" triggerIndex = "<<triggerIndex<<" size = "<<hltConfig_.size()<<std::endl;
    std::cout<<" Analyze triggerName : "<<triggerName<<std::endl;
  }
  if (triggerIndex>=n) {
    if(verbose){
      cout << "FourmuonAnalyzer4::analyzeTrigger: path "
           << triggerName << " - not found!" << endl;
    }
    return;
  }
  const unsigned int moduleIndex(hltR->index(triggerIndex));
  const unsigned int m(hltConfig_.size(triggerIndex));
  const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));
  lastTriggerModule = moduleIndex;
  if(verbose){
    cout << "FourmuonAnalyzer::analyzeTrigger: path "
         << triggerName << " [" << triggerIndex << "]" << endl;
         
    std::cout<<"  n = "<< n<<" triggerIndex = "<<triggerIndex<<" m = "<<m<<std::endl;
    std::cout<<" moduleLabels = "<<moduleLabels.size()<<" moduleIndex = "<<moduleIndex<<std::endl;
    // Results from TriggerResults product
        cout << " Trigger path status:"
         << " WasRun=" << hltR->wasrun(triggerIndex)
         << " Accept=" << hltR->accept(triggerIndex)
         << " Error =" << hltR->error(triggerIndex)
         << endl;
    cout << " Last active module - label/type: "
         << moduleLabels[moduleIndex] << "/" << hltConfig_.moduleType(moduleLabels[moduleIndex])
         << " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]"
         << endl;
  }
  assert (moduleIndex<m);
  // Results from TriggerEvent product - Attention: must look only for
  // modules actually run in this path for this event!
    std::vector < GlobalVector > passMomenta;
  for (unsigned int j=0; j<=moduleIndex; ++j) {
    const string& moduleLabel(moduleLabels[j]);
    const string  moduleType(hltConfig_.moduleType(moduleLabel));
    // check whether the module is packed up in TriggerEvent product
        const unsigned int filterIndex(hltE->filterIndex(InputTag(moduleLabel,"",hltName_)));
    if(verbose){
      std::cout<<" j = "<<j<<" modLabel/moduleType = "<<moduleLabel<<"/"<<moduleType<<" filterIndex = "<<filterIndex<<" sizeF = "<<hltE->sizeFilters()<<std::endl;
    }
    if (filterIndex<hltE->sizeFilters()) { 
          if(verbose){
        cout << " 'L3' (or 'L1', 'L2') filter in slot " << j << " - label/type " << moduleLabel << "/" << moduleType << endl;
      }
      const Vids& VIDS (hltE->filterIds(filterIndex));
      const Keys& KEYS(hltE->filterKeys(filterIndex));
      const size_type nI(VIDS.size());
      const size_type nK(KEYS.size());
      assert(nI==nK);
      const size_type n(max(nI,nK));
      if(verbose){
        cout << "   " << n  << " accepted 'L3' (or 'L1', 'L2') objects found: " << endl;
      }
      const TriggerObjectCollection& TOC(hltE->getObjects());
      for (size_type i=0; i!=n; ++i) {
        if(0==i){
          passMomenta.clear();
        }
        const TriggerObject& TO(TOC[KEYS[i]]);
        GlobalVector momentumT0(TO.px(),TO.py(),TO.pz());
        if(verbose){
          std::cout<<" i = "<<i<<" moduleLabel/moduleType : "<<moduleLabel<<"/"<<moduleType<<std::endl;
        }
        if(13==TO.id() || -13==TO.id() || 0==TO.id()){//TO.id() --> L1 Mu is always 0 (?)
          if(verbose){
            std::cout<<" current moduleType = "<<moduleType<<std::endl;
          }  
         if("HLTL1TSeed" == moduleType && "hltL1sTripleMu0orTripleMu500" == moduleLabel){ // HLT_Dimuon0_Jpsi_Muon_v5
            passMomenta.push_back(momentumT0);
            if(verbose){
               std::cout<<" L1 object found"<<std::endl;
            }
          }
          else if("HLTMuonL2FromL1TPreFilter"==moduleType){//HLT_Dimuon0_Jpsi_Muon_v5 // L2 filter names exist but never found any object
            passMomenta.push_back(momentumT0);
            if(verbose){
              std::cout<<" L2 object found"<<std::endl;
            }
          }
          else if("HLTMuonL3PreFilter"==moduleType || "HLTMuonIsoFilter"==moduleType){
            passMomenta.push_back(momentumT0);
            if(verbose){
              std::cout<<" L3 (highEff) object found"<<std::endl;
            }
          }
          else if("HLTDiMuonGlbTrkFilter"==moduleType){
            passMomenta.push_back(momentumT0);
            if(verbose){
              std::cout<<" L3 (lowEff) object found"<<std::endl;
            }
          }
          else if("HLT2MuonMuonDZ"==moduleType || "HLTMuonIsoFilter"==moduleType || ("HLTMuonL3PreFilter"==moduleType)){
            passMomenta.push_back(momentumT0);
            if(verbose){
              std::cout<<" HLT object found"<<std::endl;
            }
          }
          else if("HLTMuonDimuonL3Filter"==moduleType){//HLT_Dimuon0_Jpsi_Muon_v5
            passMomenta.push_back(momentumT0);
            if(verbose){
              std::cout<<" HLT L3 filter object found"<<std::endl;
            }
          }
          else if ("HLTDisplacedmumuFilter"==moduleType){ // HLT muons decayed from Upsilon candidate
             passMomenta.push_back(momentumT0);
             if (verbose)
               {
               std::cout<<"HLT object from (Upsilon candidate) found"<<std::endl;
                }  
             }
        }
        if(verbose){
          cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
               << TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()
               << endl;
        }
      }                    
    }
    if("HLTL1TSeed" == moduleType && "hltL1sTripleMu0orTripleMu500" == moduleLabel){ // HLT_Dimuon0_Jpsi_Muon_v5
      for(unsigned int i=0;i<passMomenta.size();++i){
        allMuL1TriggerVectors.push_back(passMomenta[i]);
      }
      if(verbose){
        std::cout<<" L1 obj FOUND; current size = "<< allMuL1TriggerVectors.size()<<std::endl;
      }
    }
    if("HLTMuonL2FromL1TPreFilter"==moduleType){//HLT_Dimuon0_Jpsi_Muon_v5
      for(unsigned int i=0;i<passMomenta.size();++i){
        allMuL2TriggerVectors.push_back(passMomenta[i]);
      }
          if(verbose){
        std::cout<<" L2 obj FOUND; current size = "<< allMuL2TriggerVectors.size()<<std::endl;
      }
    }
       if("HLTMuonL3PreFilter" == moduleType){
      for(unsigned int i=0;i<passMomenta.size();++i){
        allMuL3TriggerVectors_highEff.push_back(passMomenta[i]);
      }
          if(verbose){
        std::cout<<" L3 (highEff) obj FOUND ; current size = " <<allMuL3TriggerVectors_highEff.size()<<std::endl;
      }
    }
    if("HLTDiMuonGlbTrkFilter"==moduleType){
      for(unsigned int i=0;i<passMomenta.size();++i){
        allMuL3TriggerVectors_lowEff.push_back(passMomenta[i]);
      }
          if(verbose){
        std::cout<<" L3 (lowEff) obj FOUND ; current size = " <<allMuL3TriggerVectors_lowEff.size()<<std::endl;
      }
    }
    if("HLTMuonDimuonL3Filter" == moduleType){
      for(unsigned int i=0;i<passMomenta.size();++i){
        allMuL3TriggerVectors_highEff.push_back(passMomenta[i]);
      }
           if(verbose){
        std::cout<<" HLT L3 filter object FOUND; current size = "<<allMuL3TriggerVectors_highEff.size()<<std::endl;
      }
    }
    if("HLTMuonL3PreFilter"==moduleType){
      for(unsigned int i=0;i<passMomenta.size();++i){
        allMuHLTTriggerVectors.push_back(passMomenta[i]);
      }
            if(verbose){
        std::cout<<" HLT obj FOUND ; current size = " <<allMuHLTTriggerVectors.size()<<std::endl;
      }
    }
    if("HLTDisplacedmumuFilter"==moduleType){// HLT muons decayed from Upsilon candidate
      for(unsigned int i=0;i<passMomenta.size();++i){
         allUpsilonMuHLTTriggerVectors.push_back(passMomenta[i]);
           }
         if (verbose){
          std::cout<<"HLT object from (Upsilon candidate) found: current size ="<<allUpsilonMuHLTTriggerVectors.size()<<std::endl;
           }
        } 
    passMomenta.clear();
  }
  return;
}

bool MuMuGammaRootupler::findTrigger(edm::Handle<edm::TriggerResults> &hltR,
                          std::vector < std::string > & triggersToCheck,
                          std::vector < std::string > & triggerNameFound)
{
  triggerNameFound.clear();
  if(verbose){
    std::cout<<" findTrigger()... "<<std::endl;
    if(1==Total_events_analyzed){
      std::cout<<"   Request: "<<std::endl;
      for(unsigned int iT=0;iT<triggersToCheck.size();++iT){
        std::cout<<"       name ["<<iT<<"] = "<<triggersToCheck[iT]<<std::endl;
      }
    }
  }
  bool triggerFound = false;
  std::vector<std::string>  hlNames=triggerNames.triggerNames();
  for (uint iT=0; iT<hlNames.size(); ++iT) {
        if(verbose && 1==Total_events_analyzed){
      std::cout<<" iT = "<<iT<<" hlNames[iT] = "<<hlNames[iT]<<
        " : wasrun = "<<hltR->wasrun(iT)<<" accept = "<<
        hltR->accept(iT)<<" !error = "<<
        !hltR->error(iT)<<std::endl;
    }
   for(uint imyT = 0;imyT<triggersToCheck.size();++imyT){
      if(string::npos!=hlNames[iT].find(triggersToCheck[imyT]))
         {
        if(verbose && 1==Total_events_analyzed){
          std::cout<<" Trigger "<<hlNames[iT]<<" found to be compatible with the requested. "<<std::endl;
        }
        triggerNameFound.push_back(hlNames[iT]);
        if(triggerDecision(hltR, iT)){
          triggerFound = true;
        }
      }
    }
  }
  return triggerFound;
}

UInt_t MuMuGammaRootupler::getTriggerBits(const edm::Event& iEvent) {
	UInt_t itrigger = 0;
        // trigger collection
	edm::Handle<edm::TriggerResults> triggerResults_handle;
	iEvent.getByToken(triggerResultsTok_, triggerResults_handle);

	//   if ( triggerResults_handle.isValid() ) { 
	//     std::string testTriggerName;
	//     const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
	//     for(unsigned int trigIndex = 0; trigIndex < TheTriggerNames.size(); trigIndex++){
	//     testTriggerName = TheTriggerNames.triggerName(trigIndex);
	//     std::cout<<testTriggerName.c_str()<<std::endl;
	//     }
	//   }
	if ( triggerResults_handle.isValid() ) {
                // trigger collection 
		const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
		std::vector <unsigned int> bits_0, bits_1, bits_2, bits_3, bits_4, bits_5, bits_6, bits_7, bits_8, bits_9;
		for ( int version = 1; version<20; version ++ ) {
			std::stringstream ss0,ss1,ss2,ss3,ss4,ss5,ss6,ss7,ss8,ss9;
			ss0<<"HLT_Dimuon25_Jpsi_v"<<version;
			//ss0<<"HLT_Dimuon16_Jpsi_v"<<version;
			bits_0.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss0.str()).label().c_str()));

			ss1<<"HLT_Dimuon0_Upsilon_Muon_v"<<version;
			//ss1<<"HLT_Dimuon13_PsiPrime_v"<<version;
			bits_1.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss1.str()).label().c_str()));

			ss2<<"HLT_Dimuon13_Upsilon_v"<<version;
			bits_2.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss2.str()).label().c_str()));

			ss3<<"HLT_Dimuon0_Jpsi_Muon_v"<<version; // 2016 charmonium
			//ss3<<"HLT_Dimuon10_Jpsi_Barrel_v"<<version;
			bits_3.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss3.str()).label().c_str()));
                        ss4<<"HLT_Dimuon0_Jpsi3p5_Muon2_v"<<version; ///2017-18 charmonium
			//ss4<<"HLT_TripleMu5_v"<<version; 
			//ss4<<"HLT_Dimuon8_PsiPrime_Barrel_v"<<version;
			bits_4.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss4.str()).label().c_str()));

			ss5<<"HLT_Dimuon12_Upsilon_y1p4_v"<<version;
			//ss5<<"HLT_Dimuon8_Upsilon_Barrel_v"<<version;
			bits_5.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss5.str()).label().c_str()));

			ss6<<"HLT_Dimuon20_Jpsi_v"<<version;
			bits_6.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss6.str()).label().c_str()));

			ss7<<"HLT_Trimuon2_Upsilon5_Muon_v"<<version;
			//ss7<<"HLT_Dimuon14_Phi_Barrel_Seagulls_v"<<version;
			//ss7<<"HLT_Dimuon0_Phi_Barrel_v"<<version;
			bits_7.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss7.str()).label().c_str()));

			ss8<<"HLT_Trimuon5_3p5_2_Upsilon_Muon_v"<<version;
			//ss8<<"HLT_DoubleMu4_3_Bs_v"<<version;
			bits_8.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss8.str()).label().c_str()));

			ss9<<"HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_v"<<version;
			//ss9<<"HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v"<<version;
			bits_9.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss9.str()).label().c_str()));
		}
		for (unsigned int i=0; i<bits_0.size(); i++) {
			unsigned int bit = bits_0[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 1;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_1.size(); i++) {
			unsigned int bit = bits_1[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 2;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_2.size(); i++) {
			unsigned int bit = bits_2[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 4;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_3.size(); i++) {
			unsigned int bit = bits_3[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 8;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_4.size(); i++) {
			unsigned int bit = bits_4[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 16;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_5.size(); i++) {
			unsigned int bit = bits_5[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 32;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_6.size(); i++) {
			unsigned int bit = bits_6[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 64;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_7.size(); i++) {
			unsigned int bit = bits_7[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 128;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_8.size(); i++) {
			unsigned int bit = bits_8[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) { 
					itrigger += 256;
					break;
				}   
			}   
		}   
		for (unsigned int i=0; i<bits_9.size(); i++) {
			unsigned int bit = bits_9[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) { 
					itrigger += 512;
					break;
				}   
			}   
		} 
            } 
	return itrigger;
}

//Jinfeng 12.14
//Modify the function, the input now will be the p4 vector of the two muons in the dimuon candidate
//Since the invoking of the function is after the dimuon fit, no need to worry that it is not a dimuon candidate
bool MuMuGammaRootupler::TriggerMatch(TLorentzVector mu1p4, TLorentzVector mu2p4) {
    allTrigMuons.clear();
    if (verbose) cout << "Trigger matching for dimuon candidate" << endl;
    double reco1_eta = mu1p4.Eta();
    double reco1_phi = mu1p4.Phi();
    double reco1_pt = mu1p4.Pt();
    double reco1_mass = mu1p4.M();
    double reco2_eta = mu2p4.Eta();
    double reco2_phi = mu2p4.Phi();
    double reco2_pt = mu2p4.Pt();
    double reco2_mass = mu2p4.M();
    if (verbose) cout << "This Dimuon candidate" << " mu1pt:" << reco1_pt << " mu1eta:" << reco1_eta << " mu1phi:" << reco1_phi << endl;
    if (verbose) cout << "This Dimuon candidate" << " mu2pt:" << reco2_pt << " mu1eta:" << reco2_eta << " mu1phi:" << reco2_phi << endl;
    float dR1 = -9999.;
    float dR2 = -9999.;
    float dR1_minimum = 99;
    float dR2_minimum = 99;
    float dPt1 = 999;
    float dPt2 = 999;
    float dPt1OvP = 999;
    float dPt2OvP = 999;
    TLorentzVector TempMomenta1;
    TLorentzVector TempMomenta2;

    //Jinfeng 1.16
    float HLT_TempPt1;
    float HLT_TempPt2;

    //if (verbose) cout<<"allMuHLTTriggerVectors.size():"<<allMuHLTTriggerVectors.size()<<endl;
    if (verbose) cout << "allUpsilonMuHLTTriggerVectors.size();" << allUpsilonMuHLTTriggerVectors.size() << endl;
    //for(uint iTrig =0;iTrig<allMuHLTTriggerVectors.size();++iTrig){
    for (uint iTrig = 0; iTrig < allUpsilonMuHLTTriggerVectors.size(); ++iTrig) {
        double hlt_pt = allUpsilonMuHLTTriggerVectors[iTrig].perp();
        double hlt_eta = allUpsilonMuHLTTriggerVectors[iTrig].eta();
        double hlt_phi = allUpsilonMuHLTTriggerVectors[iTrig].phi();
        double dR1 = deltaR(reco1_eta, reco1_phi, hlt_eta, hlt_phi);
        double dR2 = deltaR(reco2_eta, reco2_phi, hlt_eta, hlt_phi);
        if (verbose) cout << "iTrig:" << iTrig << " iTrigpT:" << hlt_pt << " hlt_phi:" << hlt_phi << endl;
        if (verbose) cout << "dR1:" << dR1 << " dPt1:" << dPt1 << endl;
        if (verbose) cout << "dR2:" << dR2 << " dPt2:" << dPt2 << endl;
        if (dR1 < dR2)
        {
            if (dR1 < dR1_minimum)
            {
                dR1_minimum = dR1;
                dPt1 = std::abs(reco1_pt - hlt_pt);
                dPt1OvP = dPt1 / reco1_pt;

                if (dR1 < trg_Match_dR_cut && (dPt1 < trg_Match_dP_cut || dPt1OvP < trg_Match_dP_ratio_cut))
                {
                    if (verbose) cout << "Matching L3 sucessfull muPt1 = " << reco1_pt << " trigPt = " << hlt_pt << " dR = " << dR1 << endl;
                    TempMomenta1.SetPtEtaPhiM(reco1_pt, reco1_eta, reco1_phi, reco1_mass);
                    //1.16
                    HLT_TempPt1 = hlt_pt;
                } //Found matching to muon 1
                //Jinfeng 12.14 (comment out else line)
                //else cout << "Matching failed -> iTrig = " << hlt_pt << " eta = " << hlt_eta << " phi = " << hlt_phi << " dR =" << dR1 << endl;
            } // checking best matching object with muon 1 
        } // matching to muon1 

        if (dR2 < dR1)
        {
            if (dR2 < dR2_minimum)
            {
                dR2_minimum = dR2;
                dPt2 = std::abs(reco2_pt - hlt_pt);
                dPt2OvP = dPt2 / reco2_pt;

                if (dR2 < trg_Match_dR_cut && (dPt2 < trg_Match_dP_cut || dPt2OvP < trg_Match_dP_ratio_cut))
                {
                    if (verbose) cout << "Matching L3 sucessfull muPt2 = " << reco2_pt << " trigPt = " << hlt_pt << " dR = " << dR2 << endl;
                    TempMomenta2.SetPtEtaPhiM(reco2_pt, reco2_eta, reco2_phi, reco2_mass);
                    //1.16
                    HLT_TempPt2 = hlt_pt;
                }
                //Jinfeng 12.14 (comment out else line)
                //else cout << "Matching failed -> iTrig = " << hlt_pt << " eta = " << hlt_eta << " phi = " << hlt_phi << " dR =" << dR2 << endl;
            } //checking best matching object  with muon 2
        } // matching to muon2
    } // loop over all HLT objects

    dR1 = dR1_minimum;
    dR2 = dR2_minimum;
    v_mu1_trg_dR.push_back(dR1);
    v_mu2_trg_dR.push_back(dR2);
    if (dR1 < trg_Match_dR_cut && (dPt1 < trg_Match_dP_cut || dPt1OvP < trg_Match_dP_ratio_cut))
    {
        allTrigMuons.push_back(TempMomenta1);
        //1.16
        DuplicateHLT_Pt_1 = HLT_TempPt1;
    } // filling vector of matched muon 1
    if (dR2 < trg_Match_dR_cut && (dPt2 < trg_Match_dP_cut || dPt2OvP < trg_Match_dP_ratio_cut))
    {
        allTrigMuons.push_back(TempMomenta2);
        //1.16
        DuplicateHLT_Pt_2 = HLT_TempPt2;
    } // filling vector of matched muon 2      
    if (verbose)
    {
        cout << " AllTrigMu = " << allTrigMuons.size() << endl;
        cout << " good fourMu  run/lumi/event : " << run << "/" << lumi << "/" << event << std::endl;
        if (allTrigMuons.size() >= 2) {
            cout << "Matching good " << endl;
        }
    }
    if (allTrigMuons.size() >= 2) return true;
    else return false;
}

bool MuMuGammaRootupler::TriggerMatch_restMuons(TLorentzVector mu3p4, TLorentzVector mu4p4) {
    allRestTrigMuons.clear();
    if (verbose) cout << "Trigger matching for Rest of muons candidates" << endl;
    double reco1_eta = mu3p4.Eta();
    double reco1_phi = mu3p4.Phi();
    double reco1_pt = mu3p4.Pt();
    double reco1_mass = mu3p4.M();
    double reco2_eta = mu4p4.Eta();
    double reco2_phi = mu4p4.Phi();
    double reco2_pt = mu4p4.Pt();
    double reco2_mass = mu4p4.M();
    if (verbose) cout << "First muon candidate" << " mu1pt:" << reco1_pt << " mu1eta:" << reco1_eta << " mu1phi:" << reco1_phi << endl;
    if (verbose) cout << "Second muon candidate" << " mu2pt:" << reco2_pt << " mu2eta:" << reco2_eta << " mu2phi:" << reco2_phi << endl;
    float dR1 = -9999.;
    float dR2 = -9999.;
    float dR1_minimum = 99;
    float dR2_minimum = 99;
    float dPt1 = 999;
    float dPt2 = 999;
    float dPt1OvP = 999;
    float dPt2OvP = 999;
    TLorentzVector TempMomenta1;
    TLorentzVector TempMomenta2;
    if (verbose) cout << "allMuHLTTriggerVectors.size():" << allMuHLTTriggerVectors.size() << endl;
    for (uint iTrig = 0; iTrig < allMuHLTTriggerVectors.size(); ++iTrig) {
        double hlt_pt = allMuHLTTriggerVectors[iTrig].perp();
        double hlt_eta = allMuHLTTriggerVectors[iTrig].eta();
        double hlt_phi = allMuHLTTriggerVectors[iTrig].phi();
        double dR1 = deltaR(reco1_eta, reco1_phi, hlt_eta, hlt_phi);
        double dR2 = deltaR(reco2_eta, reco2_phi, hlt_eta, hlt_phi);

        //Jinfeng 1.16
        if ((hlt_pt == DuplicateHLT_Pt_1) || (hlt_pt == DuplicateHLT_Pt_2))
            continue;

        if (verbose) cout << "iTrig:" << iTrig << " iTrigpT:" << hlt_pt << " hlt_phi:" << hlt_phi << endl;
        if (verbose) cout << "dR1:" << dR1 << " dPt1:" << dPt1 << endl;
        if (verbose) cout << "dR2:" << dR2 << " dPt2:" << dPt2 << endl;
        if (dR1 < dR2)
        {
            if (dR1 < dR1_minimum)
            {
                dR1_minimum = dR1;
                dPt1 = std::abs(reco1_pt - hlt_pt);
                dPt1OvP = dPt1 / reco1_pt;

                if (dR1 < trg_Match_dR_cut && (dPt1 < trg_Match_dP_cut || dPt1OvP < trg_Match_dP_ratio_cut))
                {
                    if (verbose) cout << "Matching L3 sucessfull muPt1 = " << reco1_pt << " trigPt = " << hlt_pt << " dR = " << dR1 << endl;
                    TempMomenta1.SetPtEtaPhiM(reco1_pt, reco1_eta, reco1_phi, reco1_mass);
                } //Found matching to muon 1
                //Jinfeng 12.14 (comment out else line)
                //else cout << "Matching failed -> iTrig = " << hlt_pt << " eta = " << hlt_eta << " phi = " << hlt_phi << " dR =" << dR1 << endl;
            } // checking best matching object with muon 1 
        } // matching to muon1 

        if (dR2 < dR1)
        {
            if (dR2 < dR2_minimum)
            {
                dR2_minimum = dR2;
                dPt2 = std::abs(reco2_pt - hlt_pt);
                dPt2OvP = dPt2 / reco2_pt;

                if (dR2 < trg_Match_dR_cut && (dPt2 < trg_Match_dP_cut || dPt2OvP < trg_Match_dP_ratio_cut))
                {
                    if (verbose) cout << "Matching L3 sucessfull muPt2 = " << reco2_pt << " trigPt = " << hlt_pt << " dR = " << dR2 << endl;
                    TempMomenta2.SetPtEtaPhiM(reco2_pt, reco2_eta, reco2_phi, reco2_mass);
                }
                //Jinfeng 12.14 (comment out else line)
                //else cout << "Matching failed -> iTrig = " << hlt_pt << " eta = " << hlt_eta << " phi = " << hlt_phi << " dR =" << dR2 << endl;
            } //checking best matching object  with muon 2
        } // matching to muon2
    } // loop over all HLT objects

    dR1 = dR1_minimum;
    dR2 = dR2_minimum;
    fourMuFit_mu3_trg_dR.push_back(dR1);
    fourMuFit_mu4_trg_dR.push_back(dR2);
    if (dR1 < trg_Match_dR_cut && (dPt1 < trg_Match_dP_cut || dPt1OvP < trg_Match_dP_ratio_cut)) {
        allRestTrigMuons.push_back(TempMomenta1);
    } // filling vector of matched muon 1
    if (dR2 < trg_Match_dR_cut && (dPt2 < trg_Match_dP_cut || dPt2OvP < trg_Match_dP_ratio_cut)) {
        allRestTrigMuons.push_back(TempMomenta2);
    } // filling vector of matched muon 2      
    if (verbose)
    {
        cout << " All Rest Trigger matched muons size = " << allRestTrigMuons.size() << endl;
        cout << " run/lumi/event : " << run << "/" << lumi << "/" << event << std::endl;
        if (allRestTrigMuons.size() >= 1) {
            cout << "Atleast one of Rest muon Match with HLT Object " << endl;
        }
    }
    if (allRestTrigMuons.size() >= 1) return true;
    else return false;
}
// ------------ method called for each event  ------------
void MuMuGammaRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {
        using namespace trigger;
        using namespace pat;
	edm::Handle<pat::CompositeCandidateCollection> dimuons;
	iEvent.getByToken(dimuon_Label,dimuons);

	edm::Handle<pat::CompositeCandidateCollection> conversions;
	iEvent.getByToken(conversion_Label,conversions);

	edm::Handle< edm::View<pat::Muon> > muons;
	iEvent.getByToken(muon_Label, muons);

	edm::Handle<reco::VertexCollection> primaryVertices_handle;
	iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);
	reco::Vertex thePrimaryV;

	edm::Handle<reco::BeamSpot> theBeamSpot;
	iEvent.getByToken(bs_Label,theBeamSpot);
	reco::BeamSpot bs = *theBeamSpot;

	if ( primaryVertices_handle->begin() != primaryVertices_handle->end() ) {  
		thePrimaryV = reco::Vertex(*(primaryVertices_handle->begin()));
	}
	else {
		thePrimaryV = reco::Vertex(bs.position(), bs.covariance3D());
	}
	pv_x = thePrimaryV.x();
	pv_y = thePrimaryV.y();
	pv_z = thePrimaryV.z();

	/*edm::Handle<reco::ConversionCollection> conversionHandle;
	  iEvent.getByLabel(conversion_Label,conversionHandle);

	  edm::Handle<reco::PFCandidateCollection> pfcandidates;
	  iEvent.getByLabel("particleFlow",pfcandidates);
	  const reco::PFCandidateCollection pfphotons = selectPFPhotons(*pfcandidates);
	  */
         ++Total_events_analyzed;

      if(int(iEvent.id().run())!=runNumber){
      runNumber = iEvent.id().run();
      if(verbose){
        std::cout<<" New run : "<<iEvent.id().run()<<std::endl;
       }
    const edm::Run * iRun_c = &iEvent.getRun();
    edm::Run * iRun = const_cast <edm::Run*> (iRun_c);
    beginRun(*iRun, iSetup);
      }

    edm::Handle<edm::TriggerResults> hltR;
    edm::Handle<trigger::TriggerEvent> hltE;
    iEvent.getByToken(triggerEventTok_,hltE);
    iEvent.getByToken(triggerResultsTok_, hltR);
    std::vector<std::string>  hlNames;
    triggerNames = iEvent.triggerNames(*hltR);
    hlNames=triggerNames.triggerNames();
    std::vector < std::string > triggersFoundToApply;
      bool theTriggerPassed = findTrigger(hltR, triggerList, triggersFoundToApply);
      if (theTriggerPassed){
         if (verbose) cout<<" theTriggerPassed = "<<theTriggerPassed<<" checkTrigger = "<<checkTrigger<<" trigFound = "<<triggersFoundToApply.size()<<endl;
         ++Total_events_analyzed_triggered;
        }
	if (!OnlyGen_) {
		numPrimaryVertices = -1;
		if (primaryVertices_handle.isValid()) 	numPrimaryVertices = (int) primaryVertices_handle->size();

		trigger = getTriggerBits(iEvent);
		run     = iEvent.id().run();
		lumi    = iEvent.id().luminosityBlock();
		event   = iEvent.id().event();
                if (verbose) cout<<"run: "<<run<<endl;
                if (verbose) cout<<"lumi: "<<lumi<<endl;
                if (verbose) cout<<"event: "<<event<<endl;
                
}	

	//if (run < 316569)		//a temporary run number selection, 316569 is the first run of 2018A prompt reco v3
	if(true)
	{

	dimuon_pdgId = 0;
	mother_pdgId = 0;
	irank = 0;
        ups1_mass_GenMatched = -99;
        ups2_mass_GenMatched = -99;
        ups1_mass_GenMatched_ID = -99;
        ups2_mass_GenMatched_ID = -99;
        ups1_mass_GenMatched_ID_OS_VTX = -99;
        ups2_mass_GenMatched_ID_OS_VTX  = -99;
        fourmu_mass_allcuts = -99;
        fourmu_mass_trigger = -99;
        ups1_y_GenMatched = -99;
        ups2_y_GenMatched = -99;
        ups1_pt_GenMatched = -99;
        ups2_pt_GenMatched = -99;
        ups1_pt_GenMatched_ID = -99;
        ups1_y_GenMatched_ID = -99;
        ups2_pt_GenMatched_ID = -99;
        ups2_y_GenMatched_ID = -99;
        ups1_pt_GenMatched_ID_OS_VTX = -99;
        ups1_y_GenMatched_ID_OS_VTX = -99;
        ups2_pt_GenMatched_ID_OS_VTX = -99;
        ups2_y_GenMatched_ID_OS_VTX = -99;

        //Jinfeng
        temp_Lxy_noMC = -99;
        temp_Lxy_MC = -99;
        temp_CosAlpha_noMC = -99;
        temp_CosAlpha_MC = -99;
        temp_LxyPV_noMC = -99;
        temp_LxyPV_MC = -99;

        //Jinfeng 10.3
        temp_LxyPVErr_noMC = -99;
        temp_LxyPVErr_MC = -99;
        temp_LxyPVSig_noMC = -99;
        temp_LxyPVSig_MC = -99;
        temp_Distance_noMC = -99;
        temp_Distance_MC = -99;
        temp_DistanceErr_noMC = -99;
        temp_DistanceErr_MC = -99;
        temp_DistanceSig_noMC = -99;
        temp_DistanceSig_MC = -99;


        All_ups_muindex1.clear();
        All_ups_muindex2.clear();
	v_mu1Charge.clear();
	v_mu2Charge.clear();
	v_mu1_d0.clear();
	v_mu1_d0err.clear();
	v_mu2_d0.clear();
	v_mu2_d0err.clear();
	v_mu1_dz.clear();
	v_mu1_dzerr.clear();
	v_mu2_dz.clear();
	v_mu2_dzerr.clear();
	v_mu1_vz.clear();
	v_mu2_vz.clear();
        v_mu1_trg_dR.clear();
        v_mu2_trg_dR.clear();
	v_mumufit_Mass.clear();
	v_mumufit_MassErr.clear();
        v_mumufit_Mass_noMC.clear();
        v_mumufit_MassErr_noMC.clear();
        v_mumu_Pt.clear();
        v_mumu_Rapidity.clear();
	v_mumufit_VtxCL.clear();
        v_mumufit_VtxCL_noMC.clear();
	v_mumufit_VtxCL2.clear();
        v_mumufit_cTau_noMC.clear();
        v_mumufit_cTau_MC.clear();
        v_mumufit_cTauErr_noMC.clear();
        v_mumufit_cTauErr_MC.clear();

        //Jinfeng    
    v_mumufit_Lxy_noMC.clear();
    v_mumufit_Lxy_MC.clear();
    v_mumufit_CosAlpha_noMC.clear();
    v_mumufit_CosAlpha_MC.clear();
    v_mumufit_LxyPV_noMC.clear();
    v_mumufit_LxyPV_MC.clear();

    //Jinfeng 10.3
    v_mumufit_LxyPVErr_noMC.clear();
    v_mumufit_LxyPVErr_MC.clear();
    v_mumufit_LxyPVSig_noMC.clear();
    v_mumufit_LxyPVSig_MC.clear();
   
    Vector_DimuonVertex_noMC.clear();
    Vector_DimuonVertex_MC.clear();

	v_mumufit_DecayVtxX.clear();
	v_mumufit_DecayVtxY.clear();
	v_mumufit_DecayVtxZ.clear();
	v_mumufit_DecayVtxXE.clear();
	v_mumufit_DecayVtxYE.clear();
	v_mumufit_DecayVtxZE.clear();
        mu1_filtersMatched.clear();
        mu2_filtersMatched.clear();
	dimuon_p4.SetPtEtaPhiM(0,0,0,0);
	mu1_p4.SetPtEtaPhiM(0,0,0,0);
	mu2_p4.SetPtEtaPhiM(0,0,0,0);
	mu1charge = -10; 
	mu2charge = -10; 
	mu1_d0 = -10; 
	mu1_d0err = -10; 
	mu2_d0 = -10; 
	mu2_d0err = -10; 
	mu1_dz = -1000;
	mu1_dzerr = -1000;
	mu2_dz = -1000;
	mu2_dzerr = -1000;
	mu1_vz = -1000;
	mu2_vz = -1000;
	mumufit_Mass = -10; 
	mumufit_MassErr = -10; 
	mumufit_VtxCL = -10; 
	mumufit_VtxCL2 = -10; 
	mumufit_DecayVtxX = -10; 
	mumufit_DecayVtxY = -10; 
	mumufit_DecayVtxZ = -10; 
	mumufit_DecayVtxXE = -10; 
	mumufit_DecayVtxYE = -10; 
	mumufit_DecayVtxZE = -10; 
	mumufit_p4.SetPtEtaPhiM(0,0,0,0);


	fourMuFit_Mass_allComb_mix.clear();
	fourMuFit_Mass_mix = -1;
	fourMuFit_MassErr_mix = -1;
	fourMuFit_VtxX_mix = -10;
	fourMuFit_VtxY_mix = -10;
	fourMuFit_VtxZ_mix = -100;
	fourMuFit_VtxProb_mix = -1;
	fourMuFit_Chi2_mix = -10;
	fourMuFit_ndof_mix = -1;
	fourMuFit_3plus1_mix = -1;
	fourMuFit_p4_mix.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu1p4_mix.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu2p4_mix.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu3p4_mix.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu4p4_mix.SetPtEtaPhiM(0,0,0,0);
	mu3_p4_mix.SetPtEtaPhiM(0,0,0,0);
	mu4_p4_mix.SetPtEtaPhiM(0,0,0,0);
	mu3_d0_mix = -10;
	mu3_d0err_mix = -10;
	mu4_d0_mix = -10;
	mu4_d0err_mix = -10;
	mu3_dz_mix = -100;
	mu3_dzerr_mix = -100;
	mu4_dz_mix = -100;
	mu4_dzerr_mix = -100;
	mu3Charge_mix = -10;
	mu4Charge_mix = -10;
	fourMuFit_Mass_mix3evts = -1;
	fourMuFit_VtxProb_mix3evts = -1;
	fourMuFit_p4_mix3evts.SetPtEtaPhiM(0,0,0,0);
        GENfinalState=-1;
        passedFiducialSelection=false;
        GENmu_pt.clear();
        GENmu_eta.clear();
        GENmu_phi.clear();
        GENmu_mass.clear();
        GENmu_id.clear();
        GENmu_status.clear();
        GENmu_MomId.clear();
        GENmu_MomMomId.clear();
        GENups_DaughtersId.clear();
        GENups_Daughter_mupt.clear();
        GENups_Daughter_mueta.clear();
        GENups_Daughter_muphi.clear();
        GENups_Daughter_mumass.clear();
        GENups_Daughter_mustatus.clear();
        GENups_MomId.clear();
        GENups_pt.clear();
        GENups_eta.clear();
        GENups_y.clear();
        GENups_phi.clear();
        GENups_mass.clear();
        GENdimu_mass.clear();
        GENdimu_pt.clear();
        GENdimu_eta.clear();
        GENdimu_y.clear();
        GENdimu_phi.clear();
        GENX_mass.clear();
        GENX_pt.clear();
        GENX_eta.clear();
        GENX_y.clear();
        GENX_phi.clear();
        RECO_ups1_mu1_index = -99;
        RECO_ups1_mu2_index = -99;
        RECO_ups2_mu1_index = -99;
        RECO_ups2_mu2_index = -99;
	fourMuFit_Mass_allComb.clear();
	fourMuFit_Mass.clear();
	fourMuFit_MassErr.clear();
	fourMuFit_Pt.clear();
	fourMuFit_Eta.clear();
	fourMuFit_Phi.clear();
	fourMuFit_VtxX.clear();
	fourMuFit_VtxY.clear();
	fourMuFit_VtxZ.clear();
        fourMuFit_VtxXE.clear();
        fourMuFit_VtxYE.clear();
        fourMuFit_VtxZE.clear();
        fourMuFit_PVX.clear();
        fourMuFit_PVY.clear();
        fourMuFit_PVZ.clear();
        fourMuFit_PVXE.clear();
        fourMuFit_PVYE.clear();
        fourMuFit_PVZE.clear();
        fourMuFit_PV_XYE.clear();
        fourMuFit_PV_XZE.clear();
        fourMuFit_PV_YZE.clear();
        fourMuFit_PV_CL.clear();
	fourMuFit_VtxProb.clear();

        //Jinfeng_Xdistance
    fourMuFit_cosAlpha_noMC.clear();
    fourMuFit_lxy_noMC.clear();
    fourMuFit_lxyerr_noMC.clear();
    fourMuFit_lxySig_noMC.clear();
    fourMuFit_ctauPV_noMC.clear();
    fourMuFit_ctauErrPV_noMC.clear();
    fourMuFit_lxyPV_noMC.clear();
    fourMuFit_cosAlpha.clear();
    fourMuFit_lxy.clear();
    fourMuFit_lxyerr.clear();
    fourMuFit_lxySig.clear();
    fourMuFit_ctauPV.clear();
    fourMuFit_ctauErrPV.clear();
    fourMuFit_lxyPV.clear();

        fourMuFit_Double_ups_VtxProb.clear();
        fourMuFit_ups1_VtxProb.clear();
        fourMuFit_ups2_VtxProb.clear();
        fourMuFit_ups1_VtxProb_noMC.clear();
        fourMuFit_ups2_VtxProb_noMC.clear();
        fourMuFit_ups1_pt.clear();
        fourMuFit_ups2_pt.clear();
        fourMuFit_ups1_rapidity.clear();
        fourMuFit_ups2_rapidity.clear();
        fourMuFit_ups1_mass.clear();
        fourMuFit_ups2_mass.clear();
        fourMuFit_ups1_massError.clear();
        fourMuFit_ups2_massError.clear();
        fourMuFit_ups1_cTau_noMC.clear();
        fourMuFit_ups2_cTau_noMC.clear(); 
        fourMuFit_ups1_cTau_MC.clear();
        fourMuFit_ups2_cTau_MC.clear();
        fourMuFit_ups1_cTauErr_noMC.clear();
        fourMuFit_ups2_cTauErr_noMC.clear();
        fourMuFit_ups1_cTauErr_MC.clear();
        fourMuFit_ups2_cTauErr_MC.clear();

        //Jinfeng
        fourMuFit_ups1_Lxy_noMC.clear();
        fourMuFit_ups2_Lxy_noMC.clear();
        fourMuFit_ups1_Lxy_MC.clear();
        fourMuFit_ups2_Lxy_MC.clear();
        fourMuFit_ups1_CosAlpha_noMC.clear();
        fourMuFit_ups2_CosAlpha_noMC.clear();
        fourMuFit_ups1_CosAlpha_MC.clear();
        fourMuFit_ups2_CosAlpha_MC.clear();

        //Jinfeng 10.3
        fourMuFit_ups1_LxyPVErr_noMC.clear();
        fourMuFit_ups1_LxyPVErr_MC.clear();
        fourMuFit_ups2_LxyPVErr_noMC.clear();
        fourMuFit_ups2_LxyPVErr_MC.clear();
        fourMuFit_ups1_LxyPVSig_noMC.clear();
        fourMuFit_ups1_LxyPVSig_MC.clear();
        fourMuFit_ups2_LxyPVSig_noMC.clear();
        fourMuFit_ups2_LxyPVSig_MC.clear();
        fourMuFit_Distance_noMC.clear();
        fourMuFit_Distance_MC.clear();
        fourMuFit_DistanceErr_noMC.clear();
        fourMuFit_DistanceErr_MC.clear();
        fourMuFit_DistanceSig_noMC.clear();
        fourMuFit_DistanceSig_MC.clear();

        //Jinfeng 12.14
        Trigger_Matched = false;

        fourMuFit_ups1_LxyPV_noMC.clear();
        fourMuFit_ups2_LxyPV_noMC.clear();
        fourMuFit_ups1_LxyPV_MC.clear();
        fourMuFit_ups2_LxyPV_MC.clear();

        fourMuFit_ups1_VtxX.clear();
        fourMuFit_ups1_VtxY.clear();
        fourMuFit_ups1_VtxZ.clear();
        fourMuFit_ups1_VtxXE.clear();
        fourMuFit_ups1_VtxYE.clear();
        fourMuFit_ups1_VtxZE.clear();
        fourMuFit_ups2_VtxX.clear();
        fourMuFit_ups2_VtxY.clear();
        fourMuFit_ups2_VtxZ.clear();
        fourMuFit_ups2_VtxXE.clear();
        fourMuFit_ups2_VtxYE.clear();
        fourMuFit_ups2_VtxZE.clear();
        fourMuFit_wrong_ups1_mass.clear();
        fourMuFit_wrong_ups2_mass.clear();
        fourMuFit_wrong_ups1_massError.clear();
        fourMuFit_wrong_ups2_massError.clear();
	fourMuFit_ndof.clear();
	fourMuFit_mu1Pt.clear();
	fourMuFit_mu1Eta.clear();
	fourMuFit_mu1Phi.clear();
	fourMuFit_mu1E.clear();
        fourMuFit_mu1frompv.clear();
	fourMuFit_mu2Pt.clear();
	fourMuFit_mu2Eta.clear();
	fourMuFit_mu2Phi.clear();
	fourMuFit_mu2E.clear();
        fourMuFit_mu1frompv.clear();
	fourMuFit_mu3Pt.clear();
	fourMuFit_mu3Eta.clear();
	fourMuFit_mu3Phi.clear();
	fourMuFit_mu3E.clear();
        fourMuFit_mu1frompv.clear();
	fourMuFit_mu4Pt.clear();
	fourMuFit_mu4Eta.clear();
	fourMuFit_mu4Phi.clear();
	fourMuFit_mu4E.clear();
        fourMuFit_mu1frompv.clear();
        fourMuFit_mu3_trg_dR.clear();
        fourMuFit_mu4_trg_dR.clear();
        AllRecoMuons_Pt.clear();
        AllRecoMuons_Eta.clear();
        mu1_Pt.clear();
        mu1_Eta.clear();
        mu1_Phi.clear();
        mu1_E.clear();
        mu1_genindex.clear();
        mu2_Pt.clear();
        mu2_Eta.clear();
        mu2_Phi.clear();
        mu2_E.clear();
        mu2_genindex.clear();
	mu3_Pt.clear();
	mu3_Eta.clear();
	mu3_Phi.clear();
	mu3_E.clear();
        //mu3_genindex.clear();
	mu4_Pt.clear();
	mu4_Eta.clear();
	mu4_Phi.clear();
	mu4_E.clear();
        //mu4_genindex.clear();
	mu3_d0.clear();
	mu3_d0err.clear();
	mu4_d0.clear();
	mu4_d0err.clear();
	mu3_dz.clear();
	mu3_dzerr.clear();
	mu4_dz.clear();
	mu4_dzerr.clear();
	mu3Charge.clear();
	mu4Charge.clear();
	mu1_Soft.clear();
	mu2_Soft.clear();
	mu3_Soft.clear();
	mu4_Soft.clear();
        mu1_Tight.clear();
        mu2_Tight.clear();
        mu3_Tight.clear();
        mu4_Tight.clear();
	mu1_Medium.clear();
	mu2_Medium.clear();
	mu3_Medium.clear();
	mu4_Medium.clear();
	mu1_Loose.clear();
	mu2_Loose.clear();
	mu3_Loose.clear();
	mu4_Loose.clear();
	mu1_pdgID.clear();
	mu2_pdgID.clear();
	mu3_pdgID.clear();
	mu4_pdgID.clear();
        fourMuFit_mu12overlap.clear();
        fourMuFit_mu13overlap.clear();
        fourMuFit_mu14overlap.clear();
        fourMuFit_mu23overlap.clear();
        fourMuFit_mu24overlap.clear();
        fourMuFit_mu34overlap.clear();
        fourMuFit_mu12SharedSeg.clear();
        fourMuFit_mu13SharedSeg.clear();
        fourMuFit_mu14SharedSeg.clear();
        fourMuFit_mu23SharedSeg.clear();
        fourMuFit_mu24SharedSeg.clear();
        fourMuFit_mu34SharedSeg.clear();

	/*
		fourMuFit_Mass_allComb.clear();
		fourMuFit_Mass = -1;
		fourMuFit_MassErr = -1;
		fourMuFit_VtxX = -10;
		fourMuFit_VtxY = -10;
		fourMuFit_VtxZ = -100;
		fourMuFit_VtxProb = -1;
		fourMuFit_Chi2 = -10;
		fourMuFit_ndof = -1;
		fourMuFit_p4.SetPtEtaPhiM(0,0,0,0);
		fourMuFit_mu1p4.SetPtEtaPhiM(0,0,0,0);
		fourMuFit_mu2p4.SetPtEtaPhiM(0,0,0,0);
		fourMuFit_mu3p4.SetPtEtaPhiM(0,0,0,0);
		fourMuFit_mu4p4.SetPtEtaPhiM(0,0,0,0);
		mu3_p4.SetPtEtaPhiM(0,0,0,0);
		mu4_p4.SetPtEtaPhiM(0,0,0,0);
		mu3_d0 = -10;
		mu3_d0err = -10;
		mu4_d0 = -10;
		mu4_d0err = -10;
		mu3_dz = -100;
		mu3_dzerr = -100;
		mu4_dz = -100;
		mu4_dzerr = -100;
		mu3Charge = -10; 
		mu4Charge = -10; 
		mu1_Tight = -1;
		mu2_Tight = -1;
		mu3_Tight = -1;
		mu4_Tight = -1;
		mu3_pdgID = -1;
		mu4_pdgID = -1;
		theRestMuons.clear();
		upsilonMuons.clear();
		*/
	gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	gen_mu1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	gen_mu2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	dimuon_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu1_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu2_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu1Charge_bestYMass = -10; 
	mu2Charge_bestYMass = -10; 
	mu1_d0_bestYMass = -10; 
	mu1_d0err_bestYMass = -10; 
	mu2_d0_bestYMass = -10; 
	mu2_d0err_bestYMass = -10; 
	mu1_dz_bestYMass = -1000;
	mu1_dzerr_bestYMass = -1000;
	mu2_dz_bestYMass = -1000;
	mu2_dzerr_bestYMass = -1000;
	mumufit_Mass_bestYMass = -10; 
	mumufit_MassErr_bestYMass = -10; 
	mumufit_VtxCL_bestYMass = -10; 
	mumufit_VtxCL2_bestYMass = -10; 
	mumufit_DecayVtxX_bestYMass = -10; 
	mumufit_DecayVtxY_bestYMass = -10; 
	mumufit_DecayVtxZ_bestYMass = -10; 
	mumufit_DecayVtxXE_bestYMass = -10; 
	mumufit_DecayVtxYE_bestYMass = -10; 
	mumufit_DecayVtxZE_bestYMass = -10; 
	mumufit_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	bestVertex_and_bestYMass = -1;

	fourMuFit_Mass_allComb_mix_bestYMass.clear();
	fourMuFit_Mass_mix_bestYMass = -1;
	fourMuFit_MassErr_mix_bestYMass = -1;
	fourMuFit_VtxX_mix_bestYMass = -10;
	fourMuFit_VtxY_mix_bestYMass = -10;
	fourMuFit_VtxZ_mix_bestYMass = -100;
	fourMuFit_VtxProb_mix_bestYMass = -1;
	fourMuFit_Chi2_mix_bestYMass = -10;
	fourMuFit_ndof_mix_bestYMass = -1;
	fourMuFit_3plus1_mix_bestYMass = -1;
	fourMuFit_p4_mix_bestYMass.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu1p4_mix_bestYMass.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu2p4_mix_bestYMass.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu3p4_mix_bestYMass.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu4p4_mix_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu3_p4_mix_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu4_p4_mix_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu3_d0_mix_bestYMass = -10;
	mu3_d0err_mix_bestYMass = -10;
	mu4_d0_mix_bestYMass = -10;
	mu4_d0err_mix_bestYMass = -10;
	mu3_dz_mix_bestYMass = -100;
	mu3_dzerr_mix_bestYMass = -100;
	mu4_dz_mix_bestYMass = -100;
	mu4_dzerr_mix_bestYMass = -100;
	mu3Charge_mix_bestYMass = -10;
	mu4Charge_mix_bestYMass = -10;

	fourMuFit_Mass_allComb_bestYMass.clear();
	fourMuFit_Mass_bestYMass = -1;
	fourMuFit_MassErr_bestYMass = -1;
	fourMuFit_VtxX_bestYMass = -10;
	fourMuFit_VtxY_bestYMass = -10;
	fourMuFit_VtxZ_bestYMass = -100;
	fourMuFit_VtxProb_bestYMass = -1;
	fourMuFit_Chi2_bestYMass = -10;
	fourMuFit_ndof_bestYMass = -1;
	fourMuFit_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu1p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu2p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu3p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	fourMuFit_mu4p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu3_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu4_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu3_d0_bestYMass = -10;
	mu3_d0err_bestYMass = -10;
	mu4_d0_bestYMass = -10;
	mu4_d0err_bestYMass = -10;
	mu3_dz_bestYMass = -100;
	mu3_dzerr_bestYMass = -100;
	mu4_dz_bestYMass = -100;
	mu4_dzerr_bestYMass = -100;
	mu3Charge_bestYMass = -10;
	mu4Charge_bestYMass = -10;
	mu1_Tight_bestYMass = -1;
	mu2_Tight_bestYMass = -1;
	mu3_Tight_bestYMass = -1;
	mu4_Tight_bestYMass = -1;
        mu1_Medium_bestYMass = -1;
        mu2_Medium_bestYMass = -1;
        mu3_Medium_bestYMass = -1;
        mu4_Medium_bestYMass = -1;
        mu1_Loose_bestYMass = -1;
        mu2_Loose_bestYMass = -1;
        mu3_Loose_bestYMass = -1;
        mu4_Loose_bestYMass = -1;
	mu3_pdgID_bestYMass = -1;
	mu4_pdgID_bestYMass = -1;

// Trigger

  allL1TrigMuons.clear();
  allL2TrigMuons.clear();
  allL3TrigMuons.clear();
  allMuL1TriggerVectors.clear();
  allMuL2TriggerVectors.clear();
  allMuL3TriggerVectors_lowEff.clear();
  allMuL3TriggerVectors_highEff.clear();
  allMuHLTTriggerVectors.clear();
  allUpsilonMuHLTTriggerVectors.clear();

  //Jinfeng 1.16
  DuplicateHLT_Pt_1 = -1;
  DuplicateHLT_Pt_2 = -1;

  if (verbose) cout<<"triggersFoundToApply.size()"<<triggersFoundToApply.size()<<endl;
  for(unsigned int iTrig=0;iTrig<triggersFoundToApply.size();++iTrig){
    lastTriggerModule = -1;
    if (verbose) cout<<"triggersFoundToApply.at(iTrig)"<<triggersFoundToApply.at(iTrig)<<endl;
    analyzeTrigger(hltR, hltE, triggersFoundToApply.at(iTrig));
                 }
         if (verbose)cout<<"Trigger analyzed Finished"<<endl;

        edm::Handle<reco::GenParticleCollection> prunedgenParticles;
        iEvent.getByToken(prunedgenParticlesSrc_, prunedgenParticles);
        //edm::Handle<pat::PackedGenParticle> packedgenParticles;
        //iEvent.getByToken(packedgenParticlesSrc_, packedgenParticles);

	if ((isMC_ || OnlyGen_) && prunedgenParticles.isValid()) {
        if (verbose) cout<<"setting gen variables"<<endl;    
        setGENVariables(prunedgenParticles);
	}
	//if (isMC_) gen_tree->Fill();


	/*

		if ((isMC_ || OnlyGen_) && packed.isValid() && pruned.isValid()) {
		dimuon_pdgId  = 0;
		gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
		int foundit   = 0;

		for (size_t i=0; i<pruned->size(); i++) {
		int p_id = abs((*pruned)[i].pdgId());
		const reco::Candidate *aonia = &(*pruned)[i];
		if (( p_id == pdgid_ ) && (aonia->status() == 2)) {
		dimuon_pdgId = p_id;
		foundit++;
		for (size_t j=0; j<packed->size(); j++) { //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
		const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
		const reco::Candidate * d = &(*packed)[j];
		if ( motherInPrunedCollection != nullptr && (d->pdgId() == 13 ) && isAncestor(aonia , motherInPrunedCollection) ){
		gen_mu2_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
		foundit++;
		} 
		if ( motherInPrunedCollection != nullptr && (d->pdgId() == -13 ) && isAncestor(aonia , motherInPrunedCollection) ) {
		gen_mu1_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
		foundit++;
		}
		if ( foundit == 3 ) break;               
		}
		if ( foundit == 3 ) {
		gen_dimuon_p4 = gen_mu2_p4 + gen_mu1_p4;   // this should take into account FSR
		mother_pdgId  = GetAncestor(aonia)->pdgId();
		break;
		} else {
		foundit = 0;
		dimuon_pdgId = 0;
		mother_pdgId = 0;
		gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
		}            
		}  // if ( p_id
		} // for (size

	// sanity check
	//if ( ! dimuon_pdgId ) std::cout << "MuMuGammaRootupler: does not found the given decay " << run << "," << event << std::endl;
	} */ // end if isMC

	float OniaMassMax_ = OniaMassCuts_[1];
	float OniaMassMin_ = OniaMassCuts_[0];

	// Kinematic fit

	edm::ESHandle<TransientTrackBuilder> theB; 
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 
	edm::ESHandle<MagneticField> bFieldHandle;
	iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
        RefCountedKinematicTree tree_ups_part1;
        ParticleMass ups_m = upsilon_mass_;
        float ups_m_sigma = upsilon_sigma_;
        KinematicConstraint *ups_c = new MassKinematicConstraint(ups_m,ups_m_sigma);
	nGoodUpsilonCand = 0;
        nTriggeredUpsilonCand = 0;
        nTriggeredFourMuonCand = 0;

	float bestYMass = 1000;
	pat::CompositeCandidate DimuonCand_bestYMass;
        if (verbose) cout<<"dimuonCand size: "<<dimuons->size()<<endl;
        if (! OnlyGen_ && run_YY_SAF_)
        {
          //fourMuonFit(thisDimuonCand, tree_ups_part1, muons, bFieldHandle, bs, thePrimaryV);
          YY_fourMuonFit(muons, bFieldHandle, bs, thePrimaryV,primaryVertices_handle); 
         }
    

	if ( ! OnlyGen_ 
                        && dimuons.isValid() && dimuons->size() > 0 && ! run_YY_SAF_) {
		for(pat::CompositeCandidateCollection::const_iterator dimuonCand=dimuons->begin();dimuonCand!= dimuons->end(); ++dimuonCand)
		        {
			if (dimuonCand->daughter("muon1")->charge() == dimuonCand->daughter("muon2")->charge() ) continue;
			if (dimuonCand->daughter("muon1")->pt()<2.0 || dimuonCand->daughter("muon2")->pt()<2.0 ) continue;
			if (fabs(dimuonCand->daughter("muon1")->eta())>2.4|| fabs(dimuonCand->daughter("muon2")->eta())>2.4) continue;
                        if (verbose) cout<<"Muon 1 pT: "<< dimuonCand->daughter("muon1")->pt()<< "Eta: "<<dimuonCand->daughter("muon1")->eta()<<"charge: "<<dimuonCand->daughter("muon1")->charge()<<endl;
                        if (verbose) cout<<"Muon 2 pT: "<< dimuonCand->daughter("muon2")->pt()<< "Eta: "<<dimuonCand->daughter("muon2")->eta()<<"charge: "<<dimuonCand->daughter("muon2")->charge()<<endl;
			reco::TrackRef JpsiTk[2]={  //this is from Chib code
				( dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon1") ) )->innerTrack(),
				( dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon2") ) )->innerTrack()
			};

			reco::TransientTrack muon1TT(JpsiTk[0], &(*bFieldHandle));
			reco::TransientTrack muon2TT(JpsiTk[1], &(*bFieldHandle));
			reco::TrackTransientTrack muon1TTT(JpsiTk[0], &(*bFieldHandle));
			reco::TrackTransientTrack muon2TTT(JpsiTk[1], &(*bFieldHandle));

			KinematicParticleFactoryFromTransientTrack pmumuFactory;
			std::vector<RefCountedKinematicParticle> mumuParticles;
			mumuParticles.push_back(pmumuFactory.particle(muon1TT,muonMass,float(0),float(0),muonSigma));
			mumuParticles.push_back(pmumuFactory.particle(muon2TT,muonMass,float(0),float(0),muonSigma));

			KinematicParticleVertexFitter mumufitter;
			RefCountedKinematicTree mumuVertexFitTree;
                        KinematicParticleFitter csFitter;
			mumuVertexFitTree = mumufitter.fit(mumuParticles);
                        if (!(mumuVertexFitTree->isValid())) continue;
                        if (verbose) cout<<"dimuonCand have valid vertex"<<endl;
                        mumuVertexFitTree = csFitter.fit(ups_c, mumuVertexFitTree);
                        tree_ups_part1 = mumuVertexFitTree; 
                        mumuVertexFitTree->movePointerToTheTop();
			RefCountedKinematicParticle mumu_vFit_noMC;
			RefCountedKinematicVertex mumu_vFit_vertex_noMC;
                        mumu_vFit_noMC = mumuVertexFitTree->currentParticle();
                        mumu_vFit_vertex_noMC = mumuVertexFitTree->currentDecayVertex(); 
                        pat::CompositeCandidate thisDimuonCand = *dimuonCand;
                        //Jinfeng 12.14 (comment these lines out)
                        /*
                        if (isTriggerMatch_)
                        {
	                dimuon_trigger_matched = TriggerMatch(thisDimuonCand);
                        if (verbose)cout<<"Trigger matching results:"<< dimuon_trigger_matched<<endl;
	                if (!dimuon_trigger_matched) continue;
                         }
                         */
	                ++nTriggeredUpsilonCand;
                        if (verbose) cout<<"Dimuon candidate has passed the Trigger matching with dimuon mass:"<<mumu_vFit_noMC->currentState().mass()<<endl; 
                        float thisdimuon_vtxprob = ChiSquaredProbability((double)(mumu_vFit_vertex_noMC->chiSquared()),(double)(mumu_vFit_vertex_noMC->degreesOfFreedom()));    
                        if (verbose) cout<<"Dimuon mass constrained vetex fit: " <<thisdimuon_vtxprob<<endl;
			if (thisdimuon_vtxprob < upsilon_vtx_cut) continue;
			nGoodUpsilonCand++;
			fillUpsilonVector(mumuVertexFitTree,thisDimuonCand,bFieldHandle,bs);
                        if (verbose) cout<<"Filled upsilon candidate"<<nGoodUpsilonCand<<endl;
                        if (nGoodUpsilonCand==1) fillUpsilonBestVertex(mumuVertexFitTree,thisDimuonCand,bFieldHandle,bs);
                        double thisdimuon_mass = mumu_vFit_noMC->currentState().mass();
			if (best4muonCand_ == false || (best4muonCand_ == true && nGoodUpsilonCand==1)) {
				fourMuonFit(thisDimuonCand, tree_ups_part1, muons, bFieldHandle, bs, thePrimaryV);
			              }
                        //Selecting Best Candidate based on bestYMass
			if (fabs(thisdimuon_mass-upsilon_mass_) < bestYMass){
				bestYMass=fabs(thisdimuon_mass-upsilon_mass_);
                                DimuonCand_bestYMass = thisDimuonCand;
				fillUpsilonBestMass(mumuVertexFitTree,DimuonCand_bestYMass,bFieldHandle,bs);
			        }
                       //Selecting Best Candidate based on bestY vertex probablity 
                           }//atleasst one good Upsilon candidate
		         }//end of Upsilon loop
	              if (nGoodUpsilonCand>0 && best4muonCand_ == true) fourMuonFit_bestYMass(DimuonCand_bestYMass, tree_ups_part1, muons, bFieldHandle, bs, thePrimaryV);
                      if (nTriggeredUpsilonCand >= 1 && theTriggerPassed) ++Total_events_dimuon_trg_matched; //Triggered event counter
	              if (nGoodUpsilonCand>=1 && theTriggerPassed) ++Total_events_dimuon; // Total dimuon triggered events counter
                      if (fourMuFit_Mass_allComb.size()>0) ++Total_events_Fourmuon; // Total 4muon events counter
	              onia_tree->Fill();
	                  }//end run number selection
                       } //end analyze

void  MuMuGammaRootupler::fillUpsilonVector(RefCountedKinematicTree mumuVertexFitTree, pat::CompositeCandidate dimuonCand, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs) {
	mumuVertexFitTree->movePointerToTheTop();     
	RefCountedKinematicParticle mumu_vFit_noMC = mumuVertexFitTree->currentParticle();    
	RefCountedKinematicVertex mumu_vFit_vertex_noMC = mumuVertexFitTree->currentDecayVertex(); //fitted vertex is same as the commonVertex in the Onia2MuMu skim
	v_mumufit_Mass.push_back( mumu_vFit_noMC->currentState().mass());
	v_mumufit_MassErr.push_back( sqrt( mumu_vFit_noMC->currentState().kinematicParametersError().matrix()(6,6) ) ) ; 
	v_mumufit_VtxCL.push_back( ChiSquaredProbability((double)(mumu_vFit_vertex_noMC->chiSquared()),(double)(mumu_vFit_vertex_noMC->degreesOfFreedom())) );
	v_mumufit_VtxCL2.push_back( mumu_vFit_vertex_noMC->chiSquared() );
	v_mumufit_DecayVtxX.push_back( mumu_vFit_vertex_noMC->position().x() );
	v_mumufit_DecayVtxY.push_back( mumu_vFit_vertex_noMC->position().y() );
	v_mumufit_DecayVtxZ.push_back( mumu_vFit_vertex_noMC->position().z() );
	v_mumufit_DecayVtxXE.push_back( mumu_vFit_vertex_noMC->error().cxx() );
	v_mumufit_DecayVtxYE.push_back( mumu_vFit_vertex_noMC->error().cyy() );
	v_mumufit_DecayVtxZE.push_back( mumu_vFit_vertex_noMC->error().czz() );
	v_mu1Charge.push_back( dimuonCand.daughter("muon1")->charge() );
	v_mu2Charge.push_back( dimuonCand.daughter("muon2")->charge() );

	reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
	reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
	reco::TrackTransientTrack muon1TTT(muTrack1_ref, &(*bFieldHandle));
	reco::TrackTransientTrack muon2TTT(muTrack2_ref, &(*bFieldHandle));
	v_mu1_d0.push_back( -muon1TTT.dxy(bs));
	v_mu1_d0err.push_back( muon1TTT.d0Error());
	v_mu1_dz.push_back( muon1TTT.dz());
	v_mu1_dzerr.push_back( muon1TTT.dzError());
	v_mu2_d0.push_back( -muon2TTT.dxy(bs));
	v_mu2_d0err.push_back( muon2TTT.d0Error());
	v_mu2_dz.push_back( muon2TTT.dz());
	v_mu2_dzerr.push_back( muon2TTT.dzError());
}

void  MuMuGammaRootupler::fillUpsilonBestVertex(RefCountedKinematicTree mumuVertexFitTree, pat::CompositeCandidate dimuonCand, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs) {
	mumuVertexFitTree->movePointerToTheTop();     
	RefCountedKinematicParticle mumu_vFit_noMC = mumuVertexFitTree->currentParticle();    
	RefCountedKinematicVertex mumu_vFit_vertex_noMC = mumuVertexFitTree->currentDecayVertex(); //fitted vertex is same as the commonVertex in the Onia2MuMu skim
	mumufit_Mass = mumu_vFit_noMC->currentState().mass();
	mumufit_MassErr = sqrt( mumu_vFit_noMC->currentState().kinematicParametersError().matrix()(6,6) )  ;
	mumufit_VtxCL = ChiSquaredProbability((double)(mumu_vFit_vertex_noMC->chiSquared()),(double)(mumu_vFit_vertex_noMC->degreesOfFreedom())) ;
	mumufit_VtxCL2 = mumu_vFit_vertex_noMC->chiSquared() ;
	mumufit_DecayVtxX = mumu_vFit_vertex_noMC->position().x() ;
	mumufit_DecayVtxY = mumu_vFit_vertex_noMC->position().y() ;
	mumufit_DecayVtxZ = mumu_vFit_vertex_noMC->position().z() ;
	mumufit_DecayVtxXE = mumu_vFit_vertex_noMC->error().cxx() ;
	mumufit_DecayVtxYE = mumu_vFit_vertex_noMC->error().cyy() ;
	mumufit_DecayVtxZE = mumu_vFit_vertex_noMC->error().czz() ;
	mumufit_p4.SetXYZM( mumu_vFit_noMC->currentState().globalMomentum().x(), mumu_vFit_noMC->currentState().globalMomentum().y(), mumu_vFit_noMC->currentState().globalMomentum().z(), mumufit_Mass ); 

	//raw dimuon and muon
	dimuon_p4.SetPtEtaPhiM(dimuonCand.pt(), dimuonCand.eta(), dimuonCand.phi(), dimuonCand.mass());
	reco::Candidate::LorentzVector vP = dimuonCand.daughter("muon1")->p4();
	reco::Candidate::LorentzVector vM = dimuonCand.daughter("muon2")->p4();
	mu1_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
	mu2_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
	mu1charge = dimuonCand.daughter("muon1")->charge();
	mu2charge = dimuonCand.daughter("muon2")->charge();

	reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
	reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
	reco::TrackTransientTrack muon1TTT(muTrack1_ref, &(*bFieldHandle));
	reco::TrackTransientTrack muon2TTT(muTrack2_ref, &(*bFieldHandle));
	mu1_d0 = -muon1TTT.dxy(bs);
	mu1_d0err = muon1TTT.d0Error();
	mu1_dz = muon1TTT.dz();
	mu1_dzerr = muon1TTT.dzError();
	mu2_d0 = -muon2TTT.dxy(bs);
	mu2_d0err = muon2TTT.d0Error();
	mu2_dz = muon2TTT.dz();
	mu2_dzerr = muon2TTT.dzError();
}


void  MuMuGammaRootupler::fillUpsilonBestMass(RefCountedKinematicTree mumuVertexFitTree, pat::CompositeCandidate dimuonCand, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs) {
	mumuVertexFitTree->movePointerToTheTop();
	RefCountedKinematicParticle mumu_vFit_noMC = mumuVertexFitTree->currentParticle();
	RefCountedKinematicVertex mumu_vFit_vertex_noMC = mumuVertexFitTree->currentDecayVertex(); //fitted vertex is same as the commonVertex in the Onia2MuMu skim
	mumufit_Mass_bestYMass = mumu_vFit_noMC->currentState().mass();
	mumufit_MassErr_bestYMass = sqrt( mumu_vFit_noMC->currentState().kinematicParametersError().matrix()(6,6) )  ;
	mumufit_VtxCL_bestYMass = ChiSquaredProbability((double)(mumu_vFit_vertex_noMC->chiSquared()),(double)(mumu_vFit_vertex_noMC->degreesOfFreedom())) ;
	mumufit_VtxCL2_bestYMass = mumu_vFit_vertex_noMC->chiSquared() ;
	mumufit_DecayVtxX_bestYMass = mumu_vFit_vertex_noMC->position().x() ;
	mumufit_DecayVtxY_bestYMass = mumu_vFit_vertex_noMC->position().y() ;
	mumufit_DecayVtxZ_bestYMass = mumu_vFit_vertex_noMC->position().z() ;
	mumufit_DecayVtxXE_bestYMass = mumu_vFit_vertex_noMC->error().cxx() ;
	mumufit_DecayVtxYE_bestYMass = mumu_vFit_vertex_noMC->error().cyy() ;
	mumufit_DecayVtxZE_bestYMass = mumu_vFit_vertex_noMC->error().czz() ;
	mumufit_p4_bestYMass.SetXYZM( mumu_vFit_noMC->currentState().globalMomentum().x(), mumu_vFit_noMC->currentState().globalMomentum().y(), mumu_vFit_noMC->currentState().globalMomentum().z(), mumufit_Mass_bestYMass );

	//raw dimuon and muon
	dimuon_p4_bestYMass.SetPtEtaPhiM(dimuonCand.pt(), dimuonCand.eta(), dimuonCand.phi(), dimuonCand.mass());
	reco::Candidate::LorentzVector vP = dimuonCand.daughter("muon1")->p4();
	reco::Candidate::LorentzVector vM = dimuonCand.daughter("muon2")->p4();
	mu1_p4_bestYMass.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
	mu2_p4_bestYMass.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
	mu1Charge_bestYMass = dimuonCand.daughter("muon1")->charge();
	mu2Charge_bestYMass = dimuonCand.daughter("muon2")->charge();

	reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
	reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
	reco::TrackTransientTrack muon1TTT(muTrack1_ref, &(*bFieldHandle));
	reco::TrackTransientTrack muon2TTT(muTrack2_ref, &(*bFieldHandle));
	mu1_d0_bestYMass = -muon1TTT.dxy(bs);
	mu1_d0err_bestYMass = muon1TTT.d0Error();
	mu1_dz_bestYMass = muon1TTT.dz();
	mu1_dzerr_bestYMass = muon1TTT.dzError();
	mu2_d0_bestYMass = -muon2TTT.dxy(bs);
	mu2_d0err_bestYMass = muon2TTT.d0Error();
	mu2_dz_bestYMass = muon2TTT.dz();
	mu2_dzerr_bestYMass = muon2TTT.dzError();
}


// ------------ method called once each job just before starting event loop  ------------
void MuMuGammaRootupler::beginJob() {
if (verbose)
	{
	cout<<"The job Begins"<<endl;
	}
}

// ------------ method called once each job just after ending the event loop  ------------
void MuMuGammaRootupler::endJob() {

//	if (verbose)
//	{
	cout <<"Job End"<<endl;
	cout<<"Total_events_analyzed: "<<Total_events_analyzed<<endl;
	cout<<"Total_events_analyzed_triggered: "<<Total_events_analyzed_triggered<<endl;
	cout<<"Total_events_dimuon_trg_matched: "<<Total_events_dimuon_trg_matched<<endl;
        cout<<"Total_events_Fourmuon: " <<Total_events_Fourmuon<<endl;
        cout<<"Total_events_soft4mu: "<<nevent_soft4mu<<endl;
        cout<<"Total_events_nevent_4mu_pt: "<<nevent_4mu_pt<<endl;
        cout<<"Total_events_nevent_4mu_eta " <<nevent_4mu_eta<<endl;
//        }
        Total_events->SetBinContent(1,Total_events_analyzed);
			
}
// ------------ method called when ending the processing of a run  ------------
void MuMuGammaRootupler::beginRun(edm::Run & iRun, edm::EventSetup const& iSetup) {

  if(verbose){
    cout<<" New run..."<<endl;
  }
  //--- m_l1GtUtils.getL1GtRunCache(run, iSetup, true, false);
  bool hltConfigChanged;
  bool test = hltConfig_.init(iRun, iSetup, hltName_, hltConfigChanged);
  if (hltConfig_.init(iRun, iSetup, hltName_, hltConfigChanged)) {
    if(verbose){
      std::cout<<" hltConfig_.size() = "<<hltConfig_.size()<<std::endl;
    }
    // check if trigger name in (new) config
    if(verbose){
      if (triggerName_!="@") { // "@" means: analyze all triggers in config
        const unsigned int n(hltConfig_.size());
        const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
        cout<<" triggerIndex = "<<triggerIndex<<endl;
        if (triggerIndex>=n) {
          cout << "HLTEventAnalyzerAOD::beginRun:"
               << " TriggerName " << triggerName_
               << " not available in (new) config!" << endl;
          cout << "Available TriggerNames are: " << endl;
          hltConfig_.dump("Triggers");
        }
      }
      else{
        cout<<" Bad trigger name"<<endl;
      }
    }
  } else {
    cout << "beginRun:"
         << " HLT config extraction failure with process name "
         << hltName_ << endl;

  }     

}

// ------------ method called when starting to processes a run  ------------
 void MuMuGammaRootupler::endRun(edm::Run const &, edm::EventSetup const &) {
 }

// ------------ method called when starting to processes a luminosity block  ------------
void MuMuGammaRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void MuMuGammaRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuMuGammaRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

int MuMuGammaRootupler::looseMuon(edm::View<pat::Muon>::const_iterator rmu) {                                                        
	int goodLooseMuon=0;

	if(  muon::isLooseMuon(*rmu))
	goodLooseMuon = 1;
   return goodLooseMuon;
}

int MuMuGammaRootupler::mediumMuon(edm::View<pat::Muon>::const_iterator rmu) {
	int goodMediumMuon=0;

	if(  muon::isMediumMuon(*rmu))
	//bool goodGlob = rmu->isGlobalMuon() &&
	//	rmu->globalTrack()->normalizedChi2() < 3 &&
	//	rmu->combinedQuality().chi2LocalPosition < 12 &&
	//	rmu->combinedQuality().trkKink < 20;
	//if(  muon::isLooseMuon(*rmu) &&
	//		rmu->innerTrack()->validFraction() > 0.8 &&
	//		muon::segmentCompatibility(*rmu) > (goodGlob ? 0.303 : 0.451)
	//  ) 
	goodMediumMuon=1;

	return goodMediumMuon;
}

int MuMuGammaRootupler::tightMuon(edm::View<pat::Muon>::const_iterator rmu, reco::Vertex vertex) {
	int goodTightMuon=0;

	if( muon::isTightMuon(*rmu,vertex))
	//if( rmu->isGlobalMuon()
	//		&& rmu->isPFMuon()
	//		&& rmu->globalTrack()->normalizedChi2()<10.0
	//		&& rmu->globalTrack()->hitPattern().numberOfValidMuonHits()>0
	//		&& rmu->numberOfMatchedStations()>1
	//		&& fabs(rmu->muonBestTrack()->dxy( vertex.position() ))<0.2
	//		&& fabs(rmu->muonBestTrack()->dz( vertex.position() )) < 0.5
	//		&& rmu->innerTrack()->hitPattern().numberOfValidPixelHits()>0
	//		&& rmu->track()->hitPattern().trackerLayersWithMeasurement()>5
	//  ) 
	goodTightMuon = 1;

	return goodTightMuon;
}

int MuMuGammaRootupler::fourMuonMixFit(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, std::vector<pat::Muon> muons_previous1, std::vector<pat::Muon> muons_previous2, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV){
	int nGoodFourMuonMix=0;
	fourMuFit_VtxProb_mix = -1;
	for (std::vector<pat::Muon>::iterator mu3 = muons_previous1.begin(), mu3end = muons_previous1.end(); mu3 != mu3end; ++mu3){
		for (std::vector<pat::Muon>::iterator mu4 = muons_previous2.begin(), mu4end = muons_previous2.end(); mu4 != mu4end; ++mu4){
			if (mu3->charge() == mu4->charge()) continue;
			reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
			reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
			reco::TrackRef muTrack3_ref = mu3->track();
			reco::TrackRef muTrack4_ref = mu4->track();
			reco::TransientTrack muon1TT(muTrack1_ref, &(*bFieldHandle));
			reco::TransientTrack muon2TT(muTrack2_ref, &(*bFieldHandle));
			reco::TransientTrack muon3TT(muTrack3_ref, &(*bFieldHandle));
			reco::TransientTrack muon4TT(muTrack4_ref, &(*bFieldHandle));
			KinematicParticleFactoryFromTransientTrack fourMuFactory;
			std::vector<RefCountedKinematicParticle> fourMuParticles;
			fourMuParticles.push_back(fourMuFactory.particle (muon1TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon2TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));
			KinematicConstrainedVertexFitter constVertexFitter;
			//fit w/ mass constraint
			//MultiTrackKinematicConstraint *upsilon_mtc = new  TwoTrackMassKinematicConstraint(upsilon_mass_);
			//RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles,upsilon_mtc);

			//fit w/o any mass constraint
			RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles);

			if(!fourMuTree->isEmpty())
			{
				fourMuTree->movePointerToTheTop();
				RefCountedKinematicParticle fitFourMu = fourMuTree->currentParticle();
				RefCountedKinematicVertex FourMuDecayVertex = fourMuTree->currentDecayVertex();
				if (fitFourMu->currentState().isValid() &&
						ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())) > fourMuFit_VtxProb_mix)
				{ //Get chib         
					fourMuFit_Mass_mix3evts = fitFourMu->currentState().mass();
					fourMuFit_p4_mix3evts.SetXYZM(fitFourMu->currentState().kinematicParameters().momentum().x(),fitFourMu->currentState().kinematicParameters().momentum().y(),fitFourMu->currentState().kinematicParameters().momentum().z(),fourMuFit_Mass_mix3evts);
					fourMuFit_VtxProb_mix3evts = ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom()));
				}
			}
		}
	}
	return nGoodFourMuonMix;
}

int MuMuGammaRootupler::fourMuonMixFit(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, std::vector<pat::Muon> muons_previous, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV){
	int nGoodFourMuonMix=0;
	for (std::vector<pat::Muon>::iterator mu3 = muons_previous.begin(), muend = muons_previous.end(); mu3 != muend; ++mu3){
		for(std::vector<pat::Muon>::iterator mu4 = mu3+1 ; mu4 != muend; ++mu4){
			if (mu3->charge() == mu4->charge()) continue;

			TLorentzVector mu1p4,mu2p4,mu3p4,mu4p4,fourMup4;
			reco::Candidate::LorentzVector v1 = dimuonCand.daughter("muon1")->p4();
			reco::Candidate::LorentzVector v2 = dimuonCand.daughter("muon2")->p4();
			mu1p4.SetPtEtaPhiM(v1.pt(),v1.eta(),v1.phi(),v1.mass());
			mu2p4.SetPtEtaPhiM(v2.pt(),v2.eta(),v2.phi(),v2.mass());
			mu3p4.SetPtEtaPhiM(mu3->pt(), mu3->eta(), mu3->phi(), mu3->mass());
			mu4p4.SetPtEtaPhiM(mu4->pt(), mu4->eta(), mu4->phi(), mu4->mass());
			fourMup4 = mu1p4 + mu2p4 + mu3p4 + mu4p4;
			fourMuFit_Mass_allComb_mix.push_back(fourMup4.M());
		}
	}
	//std::cout<<"previousSize="<<muons_previous.size()<<", thisEventSize="<<muons->size()<<std::endl;
	int muons_previous_originalSize=muons_previous.size();

	for (edm::View<pat::Muon>::const_iterator mu = muons->begin(); mu !=  muons->end(); ++mu){
		if (mu->pt()<2.0 || fabs(mu->eta())>2.4) continue;
		if (mu-muons->begin() == dimuonCand.userInt("mu1Index"))  continue; 
		if (mu-muons->begin() == dimuonCand.userInt("mu2Index"))  continue;
		reco::GenParticleRef genMu;
		if (isMC_) genMu = mu->genParticleRef();
		if (!isMC_ || (isMC_ && !genMu.isNonnull())) {
			muons_previous.push_back(*mu);
		}
	}

	fourMuFit_VtxProb_mix = -1;
	//std::cout<<"combinedSize="<<muons_previous.size()<<std::endl;
	for (std::vector<pat::Muon>::iterator mu3 = muons_previous.begin(), muend = muons_previous.end(); mu3-muons_previous.begin() < muons_previous_originalSize; ++mu3){
		for(std::vector<pat::Muon>::iterator mu4 = mu3+1 ; mu4 != muend; ++mu4){
			if (mu3->charge() == mu4->charge()) continue;

			reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
			reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
			reco::TrackRef muTrack3_ref = mu3->track();
			reco::TrackRef muTrack4_ref = mu4->track();
			reco::TransientTrack muon1TT(muTrack1_ref, &(*bFieldHandle));
			reco::TransientTrack muon2TT(muTrack2_ref, &(*bFieldHandle));
			reco::TransientTrack muon3TT(muTrack3_ref, &(*bFieldHandle));
			reco::TransientTrack muon4TT(muTrack4_ref, &(*bFieldHandle));
			KinematicParticleFactoryFromTransientTrack fourMuFactory;
			std::vector<RefCountedKinematicParticle> fourMuParticles;
			fourMuParticles.push_back(fourMuFactory.particle (muon1TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon2TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));
			KinematicConstrainedVertexFitter constVertexFitter;
			//fit w/ mass constraint
			//MultiTrackKinematicConstraint *upsilon_mtc = new  TwoTrackMassKinematicConstraint(upsilon_mass_);
			//RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles,upsilon_mtc);

			//fit w/o any mass constraint
			RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles);

			if(!fourMuTree->isEmpty())
			{
				fourMuTree->movePointerToTheTop();
				RefCountedKinematicParticle fitFourMu = fourMuTree->currentParticle();
				RefCountedKinematicVertex FourMuDecayVertex = fourMuTree->currentDecayVertex();
				if (fitFourMu->currentState().isValid() &&
						ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())) > fourMuFit_VtxProb_mix)
				{ //Get chib         
					nGoodFourMuonMix = 1;
					fourMuFit_Mass_mix = fitFourMu->currentState().mass();
					fourMuFit_MassErr_mix = sqrt(fitFourMu->currentState().kinematicParametersError().matrix()(6,6));
					fourMuFit_p4_mix.SetXYZM(fitFourMu->currentState().kinematicParameters().momentum().x(),fitFourMu->currentState().kinematicParameters().momentum().y(),fitFourMu->currentState().kinematicParameters().momentum().z(),fourMuFit_Mass_mix);
					fourMuFit_VtxX_mix = FourMuDecayVertex->position().x();
					fourMuFit_VtxY_mix = FourMuDecayVertex->position().y();
					fourMuFit_VtxZ_mix = FourMuDecayVertex->position().z();
					fourMuFit_VtxProb_mix = ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom()));
					fourMuFit_Chi2_mix = FourMuDecayVertex->chiSquared();
					fourMuFit_ndof_mix = FourMuDecayVertex->degreesOfFreedom();
					if (mu4-muons_previous.begin() >= muons_previous_originalSize) fourMuFit_3plus1_mix = 1; 
					else fourMuFit_3plus1_mix = 0;

					//get first muon
					bool child = fourMuTree->movePointerToTheFirstChild();
					RefCountedKinematicParticle fitMu1 = fourMuTree->currentParticle();
					if(!child) break;
					double mu1M_fit_mix = fitMu1->currentState().mass();
					double mu1Px_fit_mix = fitMu1->currentState().kinematicParameters().momentum().x();
					double mu1Py_fit_mix = fitMu1->currentState().kinematicParameters().momentum().y();
					double mu1Pz_fit_mix = fitMu1->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu1p4_mix.SetXYZM( mu1Px_fit_mix, mu1Py_fit_mix, mu1Pz_fit_mix, mu1M_fit_mix );

					//get second muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu2 = fourMuTree->currentParticle();
					if(!child) break;
					double mu2M_fit_mix = fitMu2->currentState().mass();
					double mu2Px_fit_mix = fitMu2->currentState().kinematicParameters().momentum().x();
					double mu2Py_fit_mix = fitMu2->currentState().kinematicParameters().momentum().y();
					double mu2Pz_fit_mix = fitMu2->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu2p4_mix.SetXYZM( mu2Px_fit_mix, mu2Py_fit_mix, mu2Pz_fit_mix, mu2M_fit_mix );

					//get third muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu3 = fourMuTree->currentParticle();
					if(!child) break;
					double mu3M_fit_mix = fitMu3->currentState().mass();
					double mu3Px_fit_mix = fitMu3->currentState().kinematicParameters().momentum().x();
					double mu3Py_fit_mix = fitMu3->currentState().kinematicParameters().momentum().y();
					double mu3Pz_fit_mix = fitMu3->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu3p4_mix.SetXYZM( mu3Px_fit_mix, mu3Py_fit_mix, mu3Pz_fit_mix, mu3M_fit_mix );

					//get fourth muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu4 = fourMuTree->currentParticle();
					if(!child) break;
					double mu4M_fit_mix = fitMu4->currentState().mass();
					double mu4Px_fit_mix = fitMu4->currentState().kinematicParameters().momentum().x();
					double mu4Py_fit_mix = fitMu4->currentState().kinematicParameters().momentum().y();
					double mu4Pz_fit_mix = fitMu4->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu4p4_mix.SetXYZM( mu4Px_fit_mix, mu4Py_fit_mix, mu4Pz_fit_mix, mu4M_fit_mix );
					reco::Candidate::LorentzVector v3 = mu3->p4();
					reco::Candidate::LorentzVector v4 = mu4->p4();
					mu3_p4_mix.SetPtEtaPhiM(v3.pt(),v3.eta(),v3.phi(),v3.mass());
					mu4_p4_mix.SetPtEtaPhiM(v4.pt(),v4.eta(),v4.phi(),v4.mass());
					reco::TrackTransientTrack muon3TTT(muTrack3_ref, &(*bFieldHandle));
					reco::TrackTransientTrack muon4TTT(muTrack4_ref, &(*bFieldHandle));
					mu3_d0_mix = -muon3TTT.dxy(bs);
					mu3_d0err_mix = muon3TTT.d0Error();
					mu3_dz_mix = muon3TTT.dz();
					mu3_dzerr_mix = muon3TTT.dzError();
					mu4_d0_mix = -muon4TTT.dxy(bs);
					mu4_d0err_mix = muon4TTT.d0Error();
					mu4_dz_mix = muon4TTT.dz();
					mu4_dzerr_mix = muon4TTT.dzError();
					mu3Charge_mix = mu3->charge();
					mu4Charge_mix = mu4->charge();
				}
			}
		}
	}
	return nGoodFourMuonMix;
}
//void MuMuGammaRootupler::setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles, edm::Handle<pat::PackedGenParticle> packedgenParticles){
void MuMuGammaRootupler::setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles){
    reco::GenParticleCollection::const_iterator genPart;
    int j = -1;
    int nGENMuons=0;
    if (verbose) cout<<"begin looping on gen particles"<<endl;
    for(genPart = prunedgenParticles->begin(); genPart != prunedgenParticles->end(); genPart++) {
        j++;
        if (abs(genPart->pdgId())==13) {
        if (!(genPart->status()==1)) continue; // only allow status 1 gen muons
        if (!(genAna.MotherID(&prunedgenParticles->at(j))==443)) continue; // use this for standard SPS/DPS sampels for 2*Jpsi analysis
//        if (!(genAna.MotherID(&prunedgenParticles->at(j))==23)) continue; // use this line for JHU Gen samples and comment above line
        nGENMuons++;
        if (verbose) cout<<"found a gen muon: id "<<genPart->pdgId()<<" pt: "<<genPart->pt()<<" eta: "<<genPart->eta()<<" phi: "<<genPart->phi() <<endl;
           GENmu_id.push_back(genPart->pdgId() );
           GENmu_status.push_back(genPart->status());
           GENmu_pt.push_back(genPart->pt());
           GENmu_eta.push_back(genPart->eta());
           GENmu_phi.push_back(genPart->phi());
           GENmu_mass.push_back(genPart->mass());
           GENmu_MomId.push_back(genAna.MotherID(&prunedgenParticles->at(j)));
           GENmu_MomMomId.push_back(genAna.MotherMotherID(&prunedgenParticles->at(j))); 
             } // end Gen muons
           if (genPart->pdgId()==25) {
            GENX_pt.push_back(genPart->pt());
            GENX_eta.push_back(genPart->eta());
            GENX_y.push_back(genPart->rapidity());
            GENX_phi.push_back(genPart->phi());
            GENX_mass.push_back(genPart->mass());
             }
             if (verbose) cout<<"genPart id : "<<genPart->pdgId() << " status: "<<genPart->status()<<" pT:" <<genPart->pt()<< "rapidity: "<< genPart->rapidity()<<" mass:" << genPart->mass()<<endl;
//            if ((genPart->pdgId()==23) && (genPart->status()>=20 && genPart->status()<30) ) { // use this line for JHU Gen samples and comment below line
            if ((genPart->pdgId()==443) && (genPart->status()==2)) { // user this for standard SPS/DPS sampels for 2*Jpsi analysis
            const reco::Candidate *ups_dau0=genPart->daughter(0);
            int ups_dauId = fabs(ups_dau0->pdgId());
            if (verbose) cout<<"Ups daughter : "<<ups_dauId<<endl;
                int ndau = genPart->numberOfDaughters();
                if (verbose) cout<<"Ups daughter size: "<<ndau<<endl;
                for (int d=0; d<ndau; d++) {
                    const reco::Candidate *ups_dau=genPart->daughter(d);
                    if (fabs(ups_dau->pdgId()) == 13) {
                        ups_dauId = fabs(ups_dau->pdgId());
                        if (verbose) cout<<"ups Dau: "<<d<<" id: "<<fabs(ups_dau->pdgId())<<" status: "<<ups_dau->status() <<endl;
                        if (verbose) cout<<"ups Dau: "<<d<<" pt: "<<fabs(ups_dau->pt())<<" eta: "<<ups_dau->eta() << " phi: "<<ups_dau->phi()<<endl;
                        GENups_DaughtersId.push_back(fabs(ups_dau->pdgId()));
                        GENups_Daughter_mupt.push_back(ups_dau->pt());
                        GENups_Daughter_mueta.push_back(ups_dau->eta());
                        GENups_Daughter_muphi.push_back(ups_dau->phi());
                        GENups_Daughter_mumass.push_back(ups_dau->mass());
                        GENups_Daughter_mustatus.push_back(ups_dau->status());
                        GENups_Daughter_muid.push_back(ups_dau->mass());
                    }
                }
            GENups_MomId.push_back(genAna.MotherID(&prunedgenParticles->at(j)));                
            GENups_pt.push_back(genPart->pt());
            GENups_eta.push_back(genPart->eta());
            GENups_y.push_back(genPart->rapidity());
            GENups_phi.push_back(genPart->phi());
            GENups_mass.push_back(genPart->mass());
                } 
              }// finish loop on prunedgenParticles
            // pair status 1 muons to make Jpsi's
          double NGEN_muons = GENmu_pt.size();
          if (NGEN_muons < 2 ) return;
          double  minDmass=9999.0;
          double  UPS1_mass = -999.0;
          double  UPS2_mass = -999.0;
          double  UPS1_pt = -999.0;
          double  UPS2_pt = -999.0;
          double  UPS1_phi = -999.0;
          double  UPS2_phi = -999.0;
          double  UPS1_eta = -999.0;
          double  UPS2_eta = -999.0;
          double  UPS1_y = -999.0;
          double  UPS2_y = -999.0;
          GEN_ups1_mu1_index = -99;
          GEN_ups1_mu2_index = -99;
          GEN_ups2_mu1_index = -99;
          GEN_ups2_mu2_index = -99;
          for (unsigned int i=0;i<NGEN_muons;i++){
           for (unsigned int j=i+1;j<NGEN_muons;j++){
            if (GENmu_id[i]+GENmu_id[j]!=0) continue;
               TLorentzVector li, lj;
               li.SetPtEtaPhiM(GENmu_pt[i],GENmu_eta[i],GENmu_phi[i],GENmu_mass[i]);
               lj.SetPtEtaPhiM(GENmu_pt[j],GENmu_eta[j],GENmu_phi[j],GENmu_mass[j]);
               TLorentzVector lij = li+lj;
               for (unsigned int k=0;k<NGEN_muons;k++){
                 for (unsigned int l=k+1;l<NGEN_muons;l++){
                   if (i==j || i==k || i==l || j==k || j==l || l==k ) continue;
                    if (GENmu_id[k]+GENmu_id[l]!=0) continue;
                    TLorentzVector lk, ll;
                    lk.SetPtEtaPhiM(GENmu_pt[k],GENmu_eta[k],GENmu_phi[k],GENmu_mass[k]);
                    ll.SetPtEtaPhiM(GENmu_pt[l],GENmu_eta[l],GENmu_phi[l],GENmu_mass[l]); 
                    TLorentzVector lkl = lk+ll;
                    double this_massdiff = abs(lij.M()-lkl.M());
                    if (this_massdiff<minDmass)
                       {
                      minDmass=this_massdiff;
                      UPS1_mass = lij.M();
                      UPS2_mass = lkl.M();
                      UPS1_pt = lij.Pt();
                      UPS2_pt = lkl.Pt();
                      UPS1_phi = lij.Phi();
                      UPS2_phi = lkl.Phi();
                      UPS1_eta = lij.Eta();
                      UPS2_eta = lkl.Eta();
                      UPS1_y = lij.Rapidity();
                      UPS2_y = lkl.Rapidity();
                      GEN_ups1_mu1_index = i;
                      GEN_ups1_mu2_index = j;
                      GEN_ups2_mu1_index = k;
                      GEN_ups2_mu2_index = l;
                          }
                         }
                       }
                      }
                     }
            GENdimu_mass.push_back(UPS1_mass);
            GENdimu_mass.push_back(UPS2_mass);
            GENdimu_pt.push_back(UPS1_pt);
            GENdimu_pt.push_back(UPS2_pt);
            GENdimu_eta.push_back(UPS1_eta);
            GENdimu_eta.push_back(UPS2_eta);
            GENdimu_y.push_back(UPS1_y);
            GENdimu_y.push_back(UPS2_y);
            GENdimu_phi.push_back(UPS1_phi);
            GENdimu_phi.push_back(UPS2_phi);
            if (verbose) cout<<" GENdimu_mass 1 :"<<GENdimu_mass[0]<<" pt : "<<GENdimu_pt[0]<<" y: "<<GENdimu_y[0]<<endl;
            if (verbose) cout<<" GENdimu_mass 2 :"<<GENdimu_mass[1]<<" pt : "<<GENdimu_pt[1]<<" y: "<<GENdimu_y[1]<<endl;
            } // finished Gen variables 
void MuMuGammaRootupler::YY_fourMuonFit(edm::Handle< edm::View<pat::Muon> > muons, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV,  edm::Handle<reco::VertexCollection> primaryVertices_handle){
    std::vector<pat::Muon> AllMuons;
    vector<RefCountedKinematicTree> Chi;
    vector<RefCountedKinematicParticle> Chi_1;
    ParticleMass ups_m = upsilon_mass_;
    float ups_m_sigma = upsilon_sigma_;
    float omega_mass = 0.78265;
    float phi_mass = 1.019461;
    float jpsi_2s_mass = 3.686097;
    if (verbose) cout<<" Starting selecting Global or Tacker muons"<<endl;
    KinematicConstraint *ups_c = new MassKinematicConstraint(ups_m,ups_m_sigma);
    for (edm::View<pat::Muon>::const_iterator recomu = muons->begin(), muend = muons->end(); recomu != muend; ++recomu){
        if (!(recomu->isGlobalMuon() || recomu->isTrackerMuon())) continue;
        if (recomu->pt()<2.0 || fabs(recomu->eta())>2.4)  continue;
            AllMuons.push_back(*recomu);
            AllRecoMuons_Pt.push_back(recomu->pt());
            AllRecoMuons_Eta.push_back(recomu->eta());
                  } // loop all reco muons
       if (verbose) cout << "All PF muons size: " <<AllMuons.size()<<endl;
      // GEN matching of muons
      int mu_genindex[AllMuons.size()];
      float mu_gendr[AllMuons.size()]; 
      if (isMC_)
      {
       if (verbose) cout<<"begin gen matching"<<endl;
        // for each reco lepton find the nearest gen lepton with same ID
       for (unsigned int i = 0; i < AllMuons.size(); i++){
           double minDr=9999.0;
           TLorentzVector reco, gen;
           mu_genindex[i]=-1;
           mu_gendr[i]=-1;
           reco.SetPtEtaPhiM(AllMuons[i].pt(),AllMuons[i].eta(),AllMuons[i].phi(),AllMuons[i].mass());
           for (unsigned int j = 0; j < GENmu_id.size(); j++)
              {
                   gen.SetPtEtaPhiM(GENmu_pt[j],GENmu_eta[j],GENmu_phi[j],GENmu_mass[j]);
                   double thisDr = deltaR(reco.Eta(),reco.Phi(),gen.Eta(),gen.Phi());
                   if (thisDr<minDr && thisDr<0.03) {
                        mu_genindex[i]=j;
                        minDr=thisDr;
                        mu_gendr[i]=minDr;
                        }
                    } // loop over gen muons
                 if (mu_genindex[i]>-1)
                  {
                if (verbose) cout<<"Reco lepton with pT: "<<reco.Pt()<<" matched to gen lepton : "<<GENmu_pt[mu_genindex[i]]<<" deltaR: "<<minDr<<" index: "<<mu_genindex[i]<<endl;
                  }
                } // loope over all reco leptons
          } // finish GEN matching of muons
         // start Loop over all reco soft muons and make all dimuon candidates
        double N_muons = AllMuons.size();
        if (N_muons < 2 ) return;
        // if more than one reco muons matches with gen muons, consider the closest one
        // for cut flow
        int soft_temp=0;
        int temp_pt=0;
        int temp_eta=0;
        int Nmup = 0 ;
        int Nmum = 0 ;
        for (unsigned int i=0;i<N_muons;i++){
            if (AllMuons[i].charge()>0) Nmup++;
            if (AllMuons[i].charge()<0) Nmum++;
            if (muon::isSoftMuon(AllMuons[i], thePrimaryV)) soft_temp = soft_temp+1;
            if (AllMuons[i].pt()>2 ) temp_pt = temp_pt+1;
            if (abs(AllMuons[i].eta())<2.4) temp_eta = temp_eta+1;
           }
        // end cut flow 
        if (Nmup>=2 && Nmum>=2 && soft_temp==4) nevent_soft4mu++;
        if (Nmup>=2 && Nmum>=2 && soft_temp>=4 && temp_pt>=4) nevent_4mu_pt++;
        if (Nmup>=2 && Nmum>=2 && soft_temp==4 && temp_pt>=4 && temp_eta>=4) nevent_4mu_eta++;
        if (isMC_)
           {
         for (unsigned int i=0;i<N_muons;i++){
          for (unsigned int j=i+1;j<N_muons;j++){
           if ((mu_genindex[i]==mu_genindex[j]) && (mu_genindex[i]>-1&&mu_genindex[j]>-1) )
           {
           if (verbose) cout<<"Mu1 pT: "<<AllMuons[i].pt()<<" gen matched muon pT: "<<GENmu_pt[mu_genindex[i]]<<" index: "<<mu_genindex[i]<<" dr: "<<mu_gendr[i]<<endl;
           if (verbose) cout<<"Mu2 pT: "<<AllMuons[j].pt()<<" gen matched muon pT: "<<GENmu_pt[mu_genindex[j]]<<" index: "<<mu_genindex[j]<<" dr: "<<mu_gendr[j]<<endl;
           if (mu_gendr[i] < mu_gendr[j])
               {
                  mu_genindex[j] = -1;
                   }
           else {      
                  mu_genindex[i] = -1; 
                            }
                   }
                }
              }
            } // best gen reco matched muon finished
        int n_ups=0;
        for (unsigned int i=0;i<N_muons;i++)
        {
        for(unsigned int j=i+1; j<N_muons; j++){
            TLorentzVector li, lj;
            li.SetPtEtaPhiM(AllMuons[i].pt(),AllMuons[i].eta(),AllMuons[i].phi(),AllMuons[i].mass());
            lj.SetPtEtaPhiM(AllMuons[j].pt(),AllMuons[j].eta(),AllMuons[j].phi(),AllMuons[j].mass());
            TLorentzVector lilj = li+lj;
            if (verbose) cout<<"Muon 1 pT: "<<AllMuons[i].pt()<<"Muon 1 ETa: "<<AllMuons[i].eta()<<"Muon1 Charge: "<<AllMuons[i].charge()<<endl;
            if (verbose) cout<<"Muon 2 pT: "<<AllMuons[j].pt()<<"Muon 2 ETa: "<<AllMuons[j].eta()<<"Muon2 Charge: "<<AllMuons[j].charge()<<endl;
            if (verbose) cout<<"Muon 1 phi: "<<AllMuons[i].phi()<<"Muon 1 mass: "<<AllMuons[i].mass()<<endl;
            if (verbose) cout<<"Muon 2 phi: "<<AllMuons[j].phi()<<"Muon 2 mass: "<<AllMuons[j].mass()<<endl;
            if (verbose) cout<<"Dimuon mass: "<<lilj.M()<<" pT: "<<lilj.Pt()<<endl;
            if (isMC_)
             {
            if (mu_genindex[i]>-1&&mu_genindex[j]>-1)
              {
            if (verbose) cout<<"Corresping gen matched muon pT"<<GENmu_pt[mu_genindex[i]]<<" : "<<GENmu_pt[mu_genindex[j]]<<endl;
            if (verbose) cout<<"index i: "<<mu_genindex[i]<<" idnex j: "<<mu_genindex[j]<<endl;
            if ((mu_genindex[i]==GEN_ups1_mu1_index && mu_genindex[j]==GEN_ups1_mu2_index) || (mu_genindex[i]==GEN_ups1_mu2_index && mu_genindex[j]==GEN_ups1_mu1_index)) {
                ups1_mass_GenMatched = lilj.M();
                ups1_y_GenMatched = lilj.Rapidity();
                ups1_pt_GenMatched = lilj.Pt();
                RECO_ups1_mu1_index = i;
                RECO_ups1_mu2_index = j;
                   }
            if ((mu_genindex[i]==GEN_ups2_mu1_index && mu_genindex[j]==GEN_ups2_mu2_index) || (mu_genindex[i]==GEN_ups2_mu2_index && mu_genindex[j]==GEN_ups2_mu1_index)){
                ups2_mass_GenMatched = lilj.M();
                ups2_y_GenMatched = lilj.Rapidity();
                ups2_pt_GenMatched = lilj.Pt();
                RECO_ups2_mu1_index = i;
                RECO_ups2_mu2_index = j;
                  }
                } // atleast 2 reco leptons matches with gen
             } //upsilon Gen matched finished
            if ((AllMuons[i].charge()+AllMuons[j].charge())!=0) continue;
            if (AllMuons[i].pt()<2 || AllMuons[j].pt()<2) continue;
            if (abs(AllMuons[i].eta())>2.4 || abs(AllMuons[j].eta())>2.4) continue;
            if (!(muon::isSoftMuon(AllMuons[i], thePrimaryV))) continue;
            if (!(muon::isSoftMuon(AllMuons[j], thePrimaryV))) continue;
            if (isMC_)
             {
            if (mu_genindex[i]>-1&&mu_genindex[j]>-1)
              {
              if ((mu_genindex[i]==GEN_ups1_mu1_index && mu_genindex[j]==GEN_ups1_mu2_index) || (mu_genindex[i]==GEN_ups1_mu2_index && mu_genindex[j]==GEN_ups1_mu1_index)){
                ups1_mass_GenMatched_ID = lilj.M();
                ups1_pt_GenMatched_ID = lilj.Pt();
                ups1_y_GenMatched_ID = lilj.Rapidity(); 
               }
              if ((mu_genindex[i]==GEN_ups2_mu1_index && mu_genindex[j]==GEN_ups2_mu2_index) || (mu_genindex[i]==GEN_ups2_mu2_index && mu_genindex[j]==GEN_ups2_mu1_index)){  
                ups2_mass_GenMatched_ID = lilj.M();
                ups2_y_GenMatched_ID = lilj.Rapidity();
                ups2_pt_GenMatched_ID = lilj.Pt();
                }
               }// atleast 2 reco leptons matches with gen
              } //upsilon Gen matched finished
              if (verbose) cout<<" RECO_ups1_mu1_index: "<<RECO_ups1_mu1_index <<" RECO_ups1_mu2_index: "<<RECO_ups1_mu2_index<<endl;
              if (verbose) cout<<" RECO_ups2_mu1_index: "<<RECO_ups2_mu1_index <<" RECO_ups2_mu2_index: "<<RECO_ups2_mu2_index<<endl;
              if (verbose) cout<<"ups1_mass_GenMatched: "<<ups1_mass_GenMatched<<" pT: "<<ups1_pt_GenMatched<<" y: "<< ups1_y_GenMatched <<" ups2_mass_GenMatched: "<< ups2_mass_GenMatched<< " pT: "<<ups2_pt_GenMatched<<" y: "<< ups2_y_GenMatched <<endl;
              if (verbose) cout<<"ups1_mass_GenMatched_ID: "<<ups1_mass_GenMatched_ID<<" pT: "<<ups1_pt_GenMatched_ID<<" y: "<< ups1_y_GenMatched_ID<<" ups2_mass_GenMatched_ID: "<< ups2_mass_GenMatched_ID<< " pT: "<<ups2_pt_GenMatched_ID<<" y: "<< ups2_y_GenMatched_ID<<endl;
              if (verbose) cout<< "RECO OS, ID muon1 pT: " << AllMuons[i].pt() << "muon2 pT: "<<AllMuons[j].pt()<<endl;
            reco::TrackRef muTracki_ref = AllMuons[i].track();
            reco::TrackRef muTrackj_ref = AllMuons[j].track();
            reco::TransientTrack muoniTT(muTracki_ref, &(*bFieldHandle));
            reco::TransientTrack muonjTT(muTrackj_ref, &(*bFieldHandle));
            KinematicParticleFactoryFromTransientTrack mumuFactory;
            std::vector<RefCountedKinematicParticle> mumuParticles;
            mumuParticles.push_back(mumuFactory.particle (muoniTT, muonMass, float(0), float(0), muonSigma));
            mumuParticles.push_back(mumuFactory.particle (muonjTT, muonMass, float(0), float(0), muonSigma));
            KinematicParticleVertexFitter mumufitter;
            RefCountedKinematicTree mumuVertexFitTree;
            KinematicParticleFitter csFitter;
            try{mumuVertexFitTree = mumufitter.fit(mumuParticles);}
            catch (...) {
            std::cout<<"dimu fit: PerigeeKinematicState::kinematic state passed is not valid!"<<std::endl;
            continue;
            }
            if (!(mumuVertexFitTree->isValid())) continue;
            if (verbose) cout<<"This dimuon candidate has valid vertex fit."<<endl;
            // dimuon vertex fit (no mass constraint)
            RefCountedKinematicParticle mumu_vFit_noMC;
            RefCountedKinematicVertex mumu_vFit_vertex_noMC;             
            mumuVertexFitTree->movePointerToTheTop();
            mumu_vFit_noMC = mumuVertexFitTree->currentParticle();
            mumu_vFit_vertex_noMC = mumuVertexFitTree->currentDecayVertex();
            float thisdimuon_mass_noMC = mumu_vFit_noMC->currentState().mass(); 
            float thisdimuon_mass_error_noMC = sqrt( mumu_vFit_noMC->currentState().kinematicParametersError().matrix()(6,6) );
            //if (fabs(thisdimuon_mass_noMC-upsilon_mass_) > (3*1.16*sqrt(mumu_vFit_noMC->currentState().kinematicParametersError().matrix()(6,6) ))) continue;
            //if (fabs(thisdimuon_mass_noMC-upsilon_mass_) > (2.98*sqrt(mumu_vFit_noMC->currentState().kinematicParametersError().matrix()(6,6) ))) continue;
            //Jinfeng mass windows remove 2022.9.26
            if (thisdimuon_mass_noMC < 2.7 || thisdimuon_mass_noMC > 3.5) continue;
            //if (thisdimuon_mass_noMC < 2.95 || thisdimuon_mass_noMC > 3.25) continue;       
            float thisdimuon_vtxprob_noMC = ChiSquaredProbability((double)(mumu_vFit_vertex_noMC->chiSquared()),(double)(mumu_vFit_vertex_noMC->degreesOfFreedom()));
            if (verbose) cout<<" Dimuon without mass constrained vetex fit:"<<thisdimuon_vtxprob_noMC<<endl;
            if (thisdimuon_vtxprob_noMC < upsilon_vtx_cut) continue; // dimuon vtx probablity without mc 
            if (verbose) cout<<"muon candidate pass vertex Fit cut(without mass constrained)"<<endl;
            mumuVertexFitTree = csFitter.fit(ups_c, mumuVertexFitTree);
            if (!(mumuVertexFitTree->isValid())) continue;
            RefCountedKinematicParticle mumu_vFit_MC;
            RefCountedKinematicVertex mumu_vFit_vertex_MC;
            mumuVertexFitTree->movePointerToTheTop();
            mumu_vFit_MC = mumuVertexFitTree->currentParticle();
            mumu_vFit_vertex_MC = mumuVertexFitTree->currentDecayVertex();
            float thisdimuon_mass = mumu_vFit_MC->currentState().mass();
            float thisdimuon_mass_error = sqrt( mumu_vFit_MC->currentState().kinematicParametersError().matrix()(6,6) );
            if (verbose) cout << "Dimuon mass:" <<thisdimuon_mass<<endl;
            float thisdimuon_vtxprob = ChiSquaredProbability((double)(mumu_vFit_vertex_MC->chiSquared()),(double)(mumu_vFit_vertex_MC->degreesOfFreedom()));
            if (verbose) cout<<" Dimuon mass constrained vetex fit:"<<thisdimuon_vtxprob<<endl;
            //if (thisdimuon_vtxprob < upsilon_vtx_cut) continue; // dimuon vtx probablity with mc 
            if (verbose) cout<<"Dimuon candidate pass vertex Fit cut(with mass constrained)"<<endl;
            if (isMC_)
             { 
            if (mu_genindex[i]>-1&&mu_genindex[j]>-1)
              {
                if ((mu_genindex[i]==GEN_ups1_mu1_index && mu_genindex[j]==GEN_ups1_mu2_index) || (mu_genindex[i]==GEN_ups1_mu2_index && mu_genindex[j]==GEN_ups1_mu1_index))
                {
                ups1_mass_GenMatched_ID_OS_VTX = lilj.M();
                ups1_pt_GenMatched_ID_OS_VTX = lilj.Pt();
                ups1_y_GenMatched_ID_OS_VTX = lilj.Rapidity();
                 } 
                if ((mu_genindex[i]==GEN_ups2_mu1_index && mu_genindex[j]==GEN_ups2_mu2_index) || (mu_genindex[i]==GEN_ups2_mu2_index && mu_genindex[j]==GEN_ups2_mu1_index))
                {
                ups2_mass_GenMatched_ID_OS_VTX = lilj.M();
                ups2_pt_GenMatched_ID_OS_VTX = lilj.Pt();
                ups2_y_GenMatched_ID_OS_VTX = lilj.Rapidity();
                 }
                 if (verbose) cout<<"ups1_mass_GenMatched_IS_OS_vtx: "<<ups1_mass_GenMatched_ID_OS_VTX<<" pT: "<<ups1_pt_GenMatched_ID_OS_VTX<<"y: "<<ups1_y_GenMatched_ID_OS_VTX<<" ups2_mass_GenMatched_IS_OS_vtx: "<<ups2_mass_GenMatched_ID_OS_VTX<<"pT: "<<ups2_pt_GenMatched_ID_OS_VTX<<"y: "<<ups2_y_GenMatched_ID_OS_VTX<<endl;
                } // atleast 2 reco leptons matches with gen
             } //upsilon Gen matched finished 
            n_ups++;
            nGoodUpsilonCand++;
            Chi.push_back(mumuVertexFitTree);  
            All_ups_muindex1.push_back(i);
            All_ups_muindex2.push_back(j);
            mu1_genindex.push_back(mu_genindex[i]);
            mu2_genindex.push_back(mu_genindex[j]);
            //fill vectors for upsilons candidates
            v_mumufit_Mass.push_back(thisdimuon_mass);
            v_mumufit_MassErr.push_back(thisdimuon_mass_error) ;
            v_mumufit_Mass_noMC.push_back(thisdimuon_mass_noMC);
            v_mumufit_MassErr_noMC.push_back(thisdimuon_mass_error_noMC) ;
            v_mumu_Pt.push_back(lilj.Pt());
            v_mumu_Rapidity.push_back(lilj.Rapidity());       
            v_mumufit_VtxCL.push_back(thisdimuon_vtxprob);
            v_mumufit_VtxCL_noMC.push_back(thisdimuon_vtxprob_noMC);
            v_mumufit_VtxCL2.push_back( mumu_vFit_vertex_noMC->chiSquared() );
            v_mumufit_cTau_noMC.push_back((GetcTau(mumu_vFit_vertex_noMC, mumu_vFit_noMC, thePrimaryV, 0)));
            v_mumufit_cTau_MC.push_back((GetcTau(mumu_vFit_vertex_MC, mumu_vFit_MC, thePrimaryV, 1)));
            v_mumufit_cTauErr_noMC.push_back((GetcTauErr(mumu_vFit_vertex_noMC,mumu_vFit_noMC,thePrimaryV)));
            v_mumufit_cTauErr_MC.push_back((GetcTauErr(mumu_vFit_vertex_MC,mumu_vFit_MC,thePrimaryV)));


            //Jinfeng
            v_mumufit_Lxy_noMC.push_back(temp_Lxy_noMC);
            v_mumufit_Lxy_MC.push_back(temp_Lxy_MC);
            v_mumufit_CosAlpha_noMC.push_back(temp_CosAlpha_noMC);
            v_mumufit_CosAlpha_MC.push_back(temp_CosAlpha_MC);

            //Jinfeng 10.3
            v_mumufit_LxyPVErr_noMC.push_back(temp_LxyPVErr_noMC);
            v_mumufit_LxyPVErr_MC.push_back(temp_LxyPVErr_MC);
            v_mumufit_LxyPVSig_noMC.push_back(temp_LxyPVSig_noMC);
            v_mumufit_LxyPVSig_MC.push_back(temp_LxyPVSig_MC);

            Vector_DimuonVertex_noMC.push_back(mumu_vFit_vertex_noMC);
            Vector_DimuonVertex_MC.push_back(mumu_vFit_vertex_MC);

            v_mumufit_LxyPV_noMC.push_back(temp_LxyPV_noMC);
            v_mumufit_LxyPV_MC.push_back(temp_LxyPV_MC);

            v_mumufit_DecayVtxX.push_back( mumu_vFit_vertex_noMC->position().x() );
            v_mumufit_DecayVtxY.push_back( mumu_vFit_vertex_noMC->position().y() );
            v_mumufit_DecayVtxZ.push_back( mumu_vFit_vertex_noMC->position().z() );
            v_mumufit_DecayVtxXE.push_back( mumu_vFit_vertex_noMC->error().cxx() );
            v_mumufit_DecayVtxYE.push_back( mumu_vFit_vertex_noMC->error().cyy() );
            v_mumufit_DecayVtxZE.push_back( mumu_vFit_vertex_noMC->error().czz() );
            reco::TrackTransientTrack muon1TTT(muTracki_ref, &(*bFieldHandle));
            reco::TrackTransientTrack muon2TTT(muTrackj_ref, &(*bFieldHandle));
            v_mu1_d0.push_back( -muon1TTT.dxy(bs));
            v_mu1_d0err.push_back( muon1TTT.d0Error());
            v_mu1_dz.push_back( muon1TTT.dz());
            v_mu1_dzerr.push_back( muon1TTT.dzError());
            v_mu2_d0.push_back( -muon2TTT.dxy(bs));
            v_mu2_d0err.push_back( muon2TTT.d0Error());
            v_mu2_dz.push_back( muon2TTT.dz());
            v_mu2_dzerr.push_back( muon2TTT.dzError());
            v_mu1Charge.push_back( AllMuons[i].charge() );
            v_mu2Charge.push_back( AllMuons[j].charge() );
            if (verbose) cout<<" added Ups_muindex1: "<<i<<" ups_muindex2: "<<j<<endl;
            } // Loop1 over all reco muons
           } // Loop2 over all reco muons
    if (verbose) cout<<" Total dimuon candidate: "<<n_ups<<endl;
    if (n_ups < 2 ) return;
    for (int i=0; i<n_ups; i++) {
        for (int j=i+1; j<n_ups; j++) {
            int i1 = All_ups_muindex1[i]; int i2 = All_ups_muindex2[i];
            int j1 = All_ups_muindex1[j]; int j2 = All_ups_muindex2[j];
            if (verbose) cout<<"Muon 1 pT: "<<AllMuons[i1].pt()<<"Muon 1 ETa: "<<AllMuons[i1].eta()<<endl;
            if (verbose) cout<<"Muon 2 pT: "<<AllMuons[i2].pt()<<"Muon 2 ETa: "<<AllMuons[i2].eta()<<endl;
            if (verbose) cout<<"Muon 3 pT: "<<AllMuons[j1].pt()<<"Muon 3 ETa: "<<AllMuons[j1].eta()<<endl;
            if (verbose) cout<<"Muon 4 pT: "<<AllMuons[j2].pt()<<"Muon 4 ETa: "<<AllMuons[j2].eta()<<endl;
            if (i1 == j1 || i1 == j2 || i2 == j1 || i2 == j2) continue;
            TLorentzVector mu_i1, mu_i2, mu_j1, mu_j2, fourMup4;
            mu_i1.SetPtEtaPhiM(AllMuons[i1].pt(),AllMuons[i1].eta(),AllMuons[i1].phi(),AllMuons[i1].mass());
            mu_i2.SetPtEtaPhiM(AllMuons[i2].pt(),AllMuons[i2].eta(),AllMuons[i2].phi(),AllMuons[i2].mass());
            mu_j1.SetPtEtaPhiM(AllMuons[j1].pt(),AllMuons[j1].eta(),AllMuons[j1].phi(),AllMuons[j1].mass());
            mu_j2.SetPtEtaPhiM(AllMuons[j2].pt(),AllMuons[j2].eta(),AllMuons[j2].phi(),AllMuons[j2].mass());
            fourMup4 = mu_i1 + mu_i2 + mu_j1 + mu_j2;
            fourMuFit_Mass_allComb.push_back(fourMup4.M());
            reco::TrackRef muTrack1_ref = AllMuons[i1].track();
            reco::TrackRef muTrack2_ref = AllMuons[i2].track();
            reco::TrackRef muTrack3_ref = AllMuons[j1].track();
            reco::TrackRef muTrack4_ref = AllMuons[j2].track();
            reco::TransientTrack muon1TT(muTrack1_ref, &(*bFieldHandle));
            reco::TransientTrack muon2TT(muTrack2_ref, &(*bFieldHandle));
            reco::TransientTrack muon3TT(muTrack3_ref, &(*bFieldHandle));
            reco::TransientTrack muon4TT(muTrack4_ref, &(*bFieldHandle));
            KinematicParticleFactoryFromTransientTrack fourMuFactory;
            std::vector<RefCountedKinematicParticle> fourMuParticles;
            if (verbose) cout<<"after removing over laping muons"<<endl;
            if (verbose) cout<<"Muon 1 pT: "<<AllMuons[i1].pt()<<"Muon 1 ETa: "<<AllMuons[i1].eta()<<"Muon 1 charge: "<<AllMuons[i1].charge()<<endl;
            if (verbose) cout<<"Muon 2 pT: "<<AllMuons[i2].pt()<<"Muon 2 ETa: "<<AllMuons[i2].eta()<<"Muon 2 charge: "<<AllMuons[i2].charge()<<endl;
            if (verbose) cout<<"Muon 3 pT: "<<AllMuons[j1].pt()<<"Muon 3 ETa: "<<AllMuons[j1].eta()<<"Muon 3 charge: "<<AllMuons[j1].charge()<<endl;
            if (verbose) cout<<"Muon 4 pT: "<<AllMuons[j2].pt()<<"Muon 4 ETa: "<<AllMuons[j2].eta()<<"Muon 4 charge: "<<AllMuons[j2].charge()<<endl;
            fourMuParticles.push_back(fourMuFactory.particle (muon1TT, muonMass, float(0), float(0), muonSigma));
            fourMuParticles.push_back(fourMuFactory.particle (muon2TT, muonMass, float(0), float(0), muonSigma));
            fourMuParticles.push_back(fourMuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
            fourMuParticles.push_back(fourMuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));
            KinematicParticleVertexFitter  Fitter;
            RefCountedKinematicTree fourMuTree;
            try{fourMuTree = Fitter.fit(fourMuParticles);}
            catch (...) {
            std::cout<<"4mu fit: PerigeeKinematicState::kinematic state passed is not valid!"<<std::endl;
            continue;
            }
            if(fourMuTree->isEmpty()) continue;
            fourMuTree->movePointerToTheTop();
            RefCountedKinematicParticle fitFourMu = fourMuTree->currentParticle();
            RefCountedKinematicVertex FourMuDecayVertex = fourMuTree->currentDecayVertex();
            if (!(fitFourMu->currentState().isValid())) continue;
            if (verbose) cout<<"Four muon mass: "<<fitFourMu->currentState().mass()<<endl;
            if (verbose) cout<<"Four muon vertex probabilty: "<<ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom()))<<endl;
           
            //Jinfeng_Xdistance (crucial fragment)
            if (!(FourMuDecayVertex->vertexIsValid())) continue;
            TLorentzVector LV_P_noMC(fitFourMu->currentState().globalMomentum().x(), fitFourMu->currentState().globalMomentum().y(), fitFourMu->currentState().globalMomentum().z(), fitFourMu->currentState().kinematicParameters().energy());
            TVector3 vtx_noMC;
            TVector3 pvtx_noMC;
            vtx_noMC.SetXYZ(FourMuDecayVertex->position().x(), FourMuDecayVertex->position().y(), 0);
            pvtx_noMC.SetXYZ(thePrimaryV.position().x(), thePrimaryV.position().y(), 0);
            TVector3 vdiff_noMC = vtx_noMC - pvtx_noMC;
            TVector3 pperp_noMC(LV_P_noMC.Px(), LV_P_noMC.Py(), 0);
            AlgebraicVector3 vpperp_noMC;
            vpperp_noMC[0] = pperp_noMC.x();
            vpperp_noMC[1] = pperp_noMC.y();
            vpperp_noMC[2] = 0.;
            VertexDistanceXY vdistXY_noMC;
            double cosAlpha_noMC = vdiff_noMC.Dot(pperp_noMC) / (vdiff_noMC.Perp() * pperp_noMC.Perp());
            Measurement1D distXY_noMC = vdistXY_noMC.distance(Vertex(*FourMuDecayVertex), thePrimaryV);
            //double lxy_noMC = distXY_noMC.value();
            double lxy_noMC = vdiff_noMC.Perp();
            double lxyerr_noMC = distXY_noMC.error();
            double lxySig_noMC = lxy_noMC / lxyerr_noMC;
            //double ctauPV_noMC = lxy_noMC * cosAlpha_noMC * fitFourMu->currentState().mass() / pperp_noMC.Perp();
            GlobalError v1e_noMC = (Vertex(*FourMuDecayVertex)).error();
            GlobalError v2e_noMC = thePrimaryV.error();
            AlgebraicSymMatrix33 vXYe_noMC = v1e_noMC.matrix() + v2e_noMC.matrix();
            double Error_noMC = ROOT::Math::Similarity(vpperp_noMC, vXYe_noMC);
            double ctauErrPV_noMC = sqrt(Error_noMC) * fitFourMu->currentState().mass() / (pperp_noMC.Perp2());
            double lxyPV_noMC = vdiff_noMC.Dot(pperp_noMC) / pperp_noMC.Mag();
            double ctauPV_noMC = lxyPV_noMC * fitFourMu->currentState().mass() / pperp_noMC.Perp();
            fourMuFit_cosAlpha_noMC.push_back(cosAlpha_noMC);
            fourMuFit_lxy_noMC.push_back(lxy_noMC);
            fourMuFit_lxyerr_noMC.push_back(lxyerr_noMC);
            fourMuFit_lxySig_noMC.push_back(lxySig_noMC);
            fourMuFit_ctauPV_noMC.push_back(ctauPV_noMC);
            fourMuFit_ctauErrPV_noMC.push_back(ctauErrPV_noMC);
            fourMuFit_lxyPV_noMC.push_back(lxyPV_noMC);

            vector<RefCountedKinematicParticle> Chi_1;
           RefCountedKinematicTree treeupsilon_part1 = Chi[i];
           RefCountedKinematicTree treeupsilon_part2 = Chi[j];
           treeupsilon_part1->movePointerToTheTop();
           RefCountedKinematicParticle upsilon_part1 = treeupsilon_part1->currentParticle();
           Chi_1.push_back(upsilon_part1);
           treeupsilon_part2->movePointerToTheTop();
           RefCountedKinematicParticle upsilon_part2 = treeupsilon_part2->currentParticle();
           Chi_1.push_back(upsilon_part2);
           fourMuTree = Fitter.fit(Chi_1); //Two body (upsilon) Fit to get fourmu mass
           if(fourMuTree->isEmpty()) continue;
           if(!fourMuTree->isValid()) continue;
           fourMuTree->movePointerToTheTop();
           fitFourMu = fourMuTree->currentParticle();
           RefCountedKinematicVertex Double_ups_DecayVertex = fourMuTree->currentDecayVertex();
           if (!(fitFourMu->currentState().isValid())) continue;
           if (verbose) cout<<"Double ups vertex probabilty: "<<ChiSquaredProbability((double)(Double_ups_DecayVertex->chiSquared()),(double)(Double_ups_DecayVertex->degreesOfFreedom()))<<endl;
           if (verbose) cout<<"Four muon mass from two jpsi constrainted Fit: "<<fitFourMu->currentState().mass()<<endl;

           //Jinfeng 12.14
           //Jinfeng 1.16
           DuplicateHLT_Pt_1 = 0;
           DuplicateHLT_Pt_2 = 0;
           Dimuon_trigger_matched_1 = TriggerMatch(mu_i1, mu_i2);
           Rest_Muon_trigger_Matched_2 = TriggerMatch_restMuons(mu_j1, mu_j2);

           DuplicateHLT_Pt_1 = 0;
           DuplicateHLT_Pt_2 = 0;
           Dimuon_trigger_matched_2 = TriggerMatch(mu_j1, mu_j2);
           Rest_Muon_trigger_Matched_1 = TriggerMatch_restMuons(mu_i1, mu_i2);

           Trigger_Matched = ((Dimuon_trigger_matched_1 && Rest_Muon_trigger_Matched_2) || (Dimuon_trigger_matched_2 && Rest_Muon_trigger_Matched_1));

           //Jinfeng_Xdistance (crucial fragment)
           if (!(Double_ups_DecayVertex->vertexIsValid())) continue;
           TLorentzVector LV_P(fitFourMu->currentState().globalMomentum().x(), fitFourMu->currentState().globalMomentum().y(), fitFourMu->currentState().globalMomentum().z(), fitFourMu->currentState().kinematicParameters().energy());
           TVector3 vtx;
           TVector3 pvtx;
           vtx.SetXYZ(Double_ups_DecayVertex->position().x(), Double_ups_DecayVertex->position().y(), 0);
           pvtx.SetXYZ(thePrimaryV.position().x(), thePrimaryV.position().y(), 0);
           TVector3 vdiff = vtx - pvtx;
           TVector3 pperp(LV_P.Px(), LV_P.Py(), 0);
           AlgebraicVector3 vpperp;
           vpperp[0] = pperp.x();
           vpperp[1] = pperp.y();
           vpperp[2] = 0.;
           VertexDistanceXY vdistXY;
           double cosAlpha = vdiff.Dot(pperp) / (vdiff.Perp() * pperp.Perp());
           Measurement1D distXY = vdistXY.distance(Vertex(*Double_ups_DecayVertex), thePrimaryV);
           //double lxy = distXY.value();
           double lxy = vdiff.Perp();
           double lxyerr = distXY.error();
           double lxySig = lxy / lxyerr;
           //double ctauPV = lxy * cosAlpha * fitFourMu->currentState().mass() / pperp.Perp();
           GlobalError v1e = (Vertex(*Double_ups_DecayVertex)).error();
           GlobalError v2e = thePrimaryV.error();
           AlgebraicSymMatrix33 vXYe = v1e.matrix() + v2e.matrix();
           double Error = ROOT::Math::Similarity(vpperp, vXYe);
           double ctauErrPV = sqrt(Error) * fitFourMu->currentState().mass() / (pperp.Perp2());
           double lxyPV = vdiff.Dot(pperp) / pperp.Mag();
           double ctauPV = lxyPV * fitFourMu->currentState().mass() / pperp.Perp();
           fourMuFit_cosAlpha.push_back(cosAlpha);
           fourMuFit_lxy.push_back(lxy);
           fourMuFit_lxyerr.push_back(lxyerr);
           fourMuFit_lxySig.push_back(lxySig);
           fourMuFit_ctauPV.push_back(ctauPV);
           fourMuFit_ctauErrPV.push_back(ctauErrPV);
           fourMuFit_lxyPV.push_back(lxyPV);

           //Jinfeng 10.3
           GetDistance(Vector_DimuonVertex_noMC[i], Vector_DimuonVertex_noMC[j], 0);
           fourMuFit_Distance_noMC.push_back(temp_Distance_noMC);
           fourMuFit_DistanceErr_noMC.push_back(temp_DistanceErr_noMC);
           fourMuFit_DistanceSig_noMC.push_back(temp_DistanceSig_noMC);

           GetDistance(Vector_DimuonVertex_MC[i], Vector_DimuonVertex_MC[j], 1);
           fourMuFit_Distance_MC.push_back(temp_Distance_MC);
           fourMuFit_DistanceErr_MC.push_back(temp_DistanceErr_MC);
           fourMuFit_DistanceSig_MC.push_back(temp_DistanceSig_MC);
           
           fourMuFit_Mass.push_back(fitFourMu->currentState().mass());
           fourMuFit_MassErr.push_back(sqrt(fitFourMu->currentState().kinematicParametersError().matrix()(6,6)));
           TLorentzVector p4;
           p4.SetXYZM(fitFourMu->currentState().kinematicParameters().momentum().x(),fitFourMu->currentState().kinematicParameters().momentum().y(),fitFourMu->currentState().kinematicParameters().momentum().z(),fitFourMu->currentState().mass());
           // mass and mass error of wrong unselected combination, to be used to veto known resonaces 
            KinematicParticleFactoryFromTransientTrack a_MuMuFactory;
            KinematicParticleFactoryFromTransientTrack b_MuMuFactory;
            std::vector<RefCountedKinematicParticle> a_MuMuParticles;
            std::vector<RefCountedKinematicParticle> b_MuMuParticles;
            float a_wrongcomb_mass = -1.0;
            float a_wrongcomb_mass_Error = -1.0;
            float b_wrongcomb_mass = -1.0;
            float b_wrongcomb_mass_Error = -1.0;
            if (AllMuons[i1].charge()==AllMuons[j1].charge())
            {
             a_MuMuParticles.push_back(a_MuMuFactory.particle (muon1TT, muonMass, float(0), float(0), muonSigma));
             a_MuMuParticles.push_back(a_MuMuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));
             b_MuMuParticles.push_back(b_MuMuFactory.particle (muon2TT, muonMass, float(0), float(0), muonSigma));
             b_MuMuParticles.push_back(b_MuMuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
             }
            else
              {
             a_MuMuParticles.push_back(a_MuMuFactory.particle (muon1TT, muonMass, float(0), float(0), muonSigma));
             a_MuMuParticles.push_back(a_MuMuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
             b_MuMuParticles.push_back(b_MuMuFactory.particle (muon2TT, muonMass, float(0), float(0), muonSigma));
             b_MuMuParticles.push_back(b_MuMuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));
                }
            if (verbose) cout << "\tjspi1_mu1: m = "<<a_MuMuParticles[0]->currentState().mass()<<" "
            << "px = "<<a_MuMuParticles[0]->currentState().kinematicParameters().momentum().x()<<" "
            << "py = "<<a_MuMuParticles[0]->currentState().kinematicParameters().momentum().y()<<" "
            << "pz = "<<a_MuMuParticles[0]->currentState().kinematicParameters().momentum().z()<<" "
            << "x = "<<a_MuMuParticles[0]->currentState().kinematicParameters().position().x()<<" "
            << "y = "<<a_MuMuParticles[0]->currentState().kinematicParameters().position().y()<<" "
            << "z = "<<a_MuMuParticles[0]->currentState().kinematicParameters().position().z()<<" "<<endl;  
            if (verbose) cout << "\tjspi1_mu2: m = "<<a_MuMuParticles[1]->currentState().mass()<<" "
            << "px = "<<a_MuMuParticles[1]->currentState().kinematicParameters().momentum().x()<<" "
            << "py = "<<a_MuMuParticles[1]->currentState().kinematicParameters().momentum().y()<<" "
            << "pz = "<<a_MuMuParticles[1]->currentState().kinematicParameters().momentum().z()<<" "
            << "x = "<<a_MuMuParticles[1]->currentState().kinematicParameters().position().x()<<" "
            << "y = "<<a_MuMuParticles[1]->currentState().kinematicParameters().position().y()<<" "
            << "z = "<<a_MuMuParticles[1]->currentState().kinematicParameters().position().z()<<" "<<endl;
            if (verbose) cout<<"======================================"<<endl;
            if (verbose) cout << "\tjspi2_mu1: m = "<<b_MuMuParticles[0]->currentState().mass()<<" "
            << "px = "<<b_MuMuParticles[0]->currentState().kinematicParameters().momentum().x()<<" "
            << "py = "<<b_MuMuParticles[0]->currentState().kinematicParameters().momentum().y()<<" "
            << "pz = "<<b_MuMuParticles[0]->currentState().kinematicParameters().momentum().z()<<" "
            << "x = "<<b_MuMuParticles[0]->currentState().kinematicParameters().position().x()<<" "
            << "y = "<<b_MuMuParticles[0]->currentState().kinematicParameters().position().y()<<" "
            << "z = "<<b_MuMuParticles[0]->currentState().kinematicParameters().position().z()<<" "<<endl;
            if (verbose) cout << "\tjspi2_mu2: m = "<<b_MuMuParticles[1]->currentState().mass()<<" "
            << "px = "<<b_MuMuParticles[1]->currentState().kinematicParameters().momentum().x()<<" "
            << "py = "<<b_MuMuParticles[1]->currentState().kinematicParameters().momentum().y()<<" "
            << "pz = "<<b_MuMuParticles[1]->currentState().kinematicParameters().momentum().z()<<" "
            << "x = "<<b_MuMuParticles[1]->currentState().kinematicParameters().position().x()<<" "
            << "y = "<<b_MuMuParticles[1]->currentState().kinematicParameters().position().y()<<" "
            << "z = "<<b_MuMuParticles[1]->currentState().kinematicParameters().position().z()<<" "<<endl;
            KinematicParticleVertexFitter  a_Fitter;
            RefCountedKinematicTree a_MuMuTree = a_Fitter.fit(a_MuMuParticles);
            if (!(a_MuMuTree->isEmpty()))
            {
             a_MuMuTree->movePointerToTheTop();
             RefCountedKinematicParticle a_MuMu = a_MuMuTree->currentParticle();
             a_wrongcomb_mass = a_MuMu->currentState().mass();
             a_wrongcomb_mass_Error = sqrt(a_MuMu->currentState().kinematicParametersError().matrix()(6,6));
             }
            KinematicParticleVertexFitter  b_Fitter;
            RefCountedKinematicTree b_MuMuTree = b_Fitter.fit(b_MuMuParticles);
            if (!(b_MuMuTree->isEmpty()))
            {
             b_MuMuTree->movePointerToTheTop();
             RefCountedKinematicParticle b_MuMu = b_MuMuTree->currentParticle();
             b_wrongcomb_mass = b_MuMu->currentState().mass();
             b_wrongcomb_mass_Error = sqrt(b_MuMu->currentState().kinematicParametersError().matrix()(6,6));
             }           
           // finished, mass and mass error of wrong unselected combination
           fourMuFit_Pt.push_back(p4.Pt());
           fourMuFit_Eta.push_back(p4.Eta());
           fourMuFit_Phi.push_back(p4.Phi());
           fourMuFit_VtxX.push_back(FourMuDecayVertex->position().x());
           fourMuFit_VtxY.push_back(FourMuDecayVertex->position().y());
           fourMuFit_VtxZ.push_back(FourMuDecayVertex->position().z());
           fourMuFit_VtxXE.push_back(FourMuDecayVertex->error().cxx());
           fourMuFit_VtxYE.push_back(FourMuDecayVertex->error().cyy());
           fourMuFit_VtxZE.push_back(FourMuDecayVertex->error().czz());
           float temp_fourMuFit_VtxProb = ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom()));
           fourMuFit_VtxProb.push_back(temp_fourMuFit_VtxProb);
           fourMuFit_Double_ups_VtxProb.push_back(ChiSquaredProbability((double)(Double_ups_DecayVertex->chiSquared()),(double)(Double_ups_DecayVertex->degreesOfFreedom())));
           fourMuFit_ups1_VtxProb.push_back(v_mumufit_VtxCL[i]);
           fourMuFit_ups2_VtxProb.push_back(v_mumufit_VtxCL[j]);
           fourMuFit_ups1_VtxProb_noMC.push_back(v_mumufit_VtxCL_noMC[i]);
           fourMuFit_ups2_VtxProb_noMC.push_back(v_mumufit_VtxCL_noMC[j]);
           fourMuFit_ups1_pt.push_back(v_mumu_Pt[i]);
           fourMuFit_ups2_pt.push_back(v_mumu_Pt[j]);
           fourMuFit_ups1_rapidity.push_back(v_mumu_Rapidity[i]);
           fourMuFit_ups2_rapidity.push_back(v_mumu_Rapidity[j]);
           fourMuFit_ups1_mass.push_back(v_mumufit_Mass_noMC[i]);
           fourMuFit_ups2_mass.push_back(v_mumufit_Mass_noMC[j]);
           fourMuFit_ups1_massError.push_back(v_mumufit_MassErr_noMC[i]);
           fourMuFit_ups2_massError.push_back(v_mumufit_MassErr_noMC[j]);
           fourMuFit_ups1_cTau_noMC.push_back(v_mumufit_cTau_noMC[i]);
           fourMuFit_ups2_cTau_noMC.push_back(v_mumufit_cTau_noMC[j]);
           fourMuFit_ups1_cTau_MC.push_back(v_mumufit_cTau_MC[i]);
           fourMuFit_ups2_cTau_MC.push_back(v_mumufit_cTau_MC[j]);
           fourMuFit_ups1_cTauErr_noMC.push_back(v_mumufit_cTauErr_noMC[i]);
           fourMuFit_ups2_cTauErr_noMC.push_back(v_mumufit_cTauErr_noMC[j]);
           fourMuFit_ups1_cTauErr_MC.push_back(v_mumufit_cTauErr_MC[i]);
           fourMuFit_ups2_cTauErr_MC.push_back(v_mumufit_cTauErr_MC[j]);

           //Jinfeng
           fourMuFit_ups1_Lxy_noMC.push_back(v_mumufit_Lxy_noMC[i]);
           fourMuFit_ups2_Lxy_noMC.push_back(v_mumufit_Lxy_noMC[j]);
           fourMuFit_ups1_Lxy_MC.push_back(v_mumufit_Lxy_MC[i]);
           fourMuFit_ups2_Lxy_MC.push_back(v_mumufit_Lxy_MC[j]);
           fourMuFit_ups1_CosAlpha_noMC.push_back(v_mumufit_CosAlpha_noMC[i]);
           fourMuFit_ups2_CosAlpha_noMC.push_back(v_mumufit_CosAlpha_noMC[j]);
           fourMuFit_ups1_CosAlpha_MC.push_back(v_mumufit_CosAlpha_MC[i]);
           fourMuFit_ups2_CosAlpha_MC.push_back(v_mumufit_CosAlpha_MC[j]);

           //Jinfeng 10.3
           fourMuFit_ups1_LxyPVErr_noMC.push_back(v_mumufit_LxyPVErr_noMC[i]);
           fourMuFit_ups2_LxyPVErr_noMC.push_back(v_mumufit_LxyPVErr_noMC[j]);
           fourMuFit_ups1_LxyPVErr_MC.push_back(v_mumufit_LxyPVErr_MC[i]);
           fourMuFit_ups2_LxyPVErr_MC.push_back(v_mumufit_LxyPVErr_MC[j]);
           fourMuFit_ups1_LxyPVSig_noMC.push_back(v_mumufit_LxyPVSig_noMC[i]);
           fourMuFit_ups2_LxyPVSig_noMC.push_back(v_mumufit_LxyPVSig_noMC[j]);
           fourMuFit_ups1_LxyPVSig_MC.push_back(v_mumufit_LxyPVSig_MC[i]);
           fourMuFit_ups2_LxyPVSig_MC.push_back(v_mumufit_LxyPVSig_MC[j]);

           fourMuFit_ups1_LxyPV_noMC.push_back(v_mumufit_LxyPV_noMC[i]);
           fourMuFit_ups2_LxyPV_noMC.push_back(v_mumufit_LxyPV_noMC[j]);
           fourMuFit_ups1_LxyPV_MC.push_back(v_mumufit_LxyPV_MC[i]);
           fourMuFit_ups2_LxyPV_MC.push_back(v_mumufit_LxyPV_MC[j]);

           fourMuFit_ups1_VtxX.push_back(v_mumufit_DecayVtxX[i]);
           fourMuFit_ups1_VtxY.push_back(v_mumufit_DecayVtxY[i]);
           fourMuFit_ups1_VtxZ.push_back(v_mumufit_DecayVtxZ[i]);
           fourMuFit_ups1_VtxXE.push_back(v_mumufit_DecayVtxXE[i]);
           fourMuFit_ups1_VtxYE.push_back(v_mumufit_DecayVtxYE[i]);
           fourMuFit_ups1_VtxZE.push_back(v_mumufit_DecayVtxZE[i]);
           fourMuFit_ups2_VtxX.push_back(v_mumufit_DecayVtxX[j]);
           fourMuFit_ups2_VtxY.push_back(v_mumufit_DecayVtxY[j]);
           fourMuFit_ups2_VtxZ.push_back(v_mumufit_DecayVtxZ[j]);
           fourMuFit_ups2_VtxXE.push_back(v_mumufit_DecayVtxXE[j]);
           fourMuFit_ups2_VtxYE.push_back(v_mumufit_DecayVtxYE[j]);
           fourMuFit_ups2_VtxZE.push_back(v_mumufit_DecayVtxZE[j]);
           fourMuFit_ndof.push_back(FourMuDecayVertex->degreesOfFreedom());
           fourMuFit_Double_ups_Chi2.push_back(Double_ups_DecayVertex->chiSquared());
           fourMuFit_Double_ups_ndof.push_back(Double_ups_DecayVertex->degreesOfFreedom());
           // wrong combination
           if (verbose) cout<<"a_wrongcomb_mass: "<<a_wrongcomb_mass<<endl;
           if (verbose) cout<<"a_wrongcomb_mass_Error: "<<a_wrongcomb_mass_Error<<endl;
           if (verbose) cout<<"b_wrongcomb_mass: "<<b_wrongcomb_mass<<endl;
           if (verbose) cout<<"b_wrongcomb_mass_Error: "<<b_wrongcomb_mass_Error<<endl;
           fourMuFit_wrong_ups1_mass.push_back(a_wrongcomb_mass);
           fourMuFit_wrong_ups2_mass.push_back(b_wrongcomb_mass);
           fourMuFit_wrong_ups1_massError.push_back(a_wrongcomb_mass_Error);
           fourMuFit_wrong_ups2_massError.push_back(b_wrongcomb_mass_Error);
           //get first muon //Childs are taken from JPsi 1 and 2
           bool child = treeupsilon_part1->movePointerToTheFirstChild();
           RefCountedKinematicParticle fitMu1 = treeupsilon_part1->currentParticle();
           if(!child) break;
           TLorentzVector mu1p4;
           mu1p4.SetXYZM( fitMu1->currentState().kinematicParameters().momentum().x(), fitMu1->currentState().kinematicParameters().momentum().y(),fitMu1->currentState().kinematicParameters().momentum().z(), fitMu1->currentState().mass() );
           fourMuFit_mu1Pt.push_back(mu1p4.Pt());
           fourMuFit_mu1Eta.push_back(mu1p4.Eta());
           fourMuFit_mu1Phi.push_back(mu1p4.Phi());
           fourMuFit_mu1E.push_back(mu1p4.E());
           if (verbose) cout<<"post Fit mu1pt: "<<mu1p4.Pt()<<"mu1eta: "<<mu1p4.Eta()<<"mu1phi: "<<mu1p4.Phi()<<"mu1E: "<<mu1p4.E()<<endl;
           //get second muon
           child = treeupsilon_part1->movePointerToTheNextChild();
           RefCountedKinematicParticle fitMu2 = treeupsilon_part1->currentParticle();
           if(!child) break;
           TLorentzVector mu2p4;
           mu2p4.SetXYZM( fitMu2->currentState().kinematicParameters().momentum().x(), fitMu2->currentState().kinematicParameters().momentum().y(),fitMu2->currentState().kinematicParameters().momentum().z(), fitMu2->currentState().mass() );
           fourMuFit_mu2Pt.push_back(mu2p4.Pt());
           fourMuFit_mu2Eta.push_back(mu2p4.Eta());
           fourMuFit_mu2Phi.push_back(mu2p4.Phi());
           fourMuFit_mu2E.push_back(mu2p4.E());
           if (verbose) cout<<"post Fit mu2pt: "<<mu2p4.Pt()<<"mu2eta: "<<mu2p4.Eta()<<"mu2phi: "<<mu2p4.Phi()<<"mu2E: "<<mu2p4.E()<<endl;
           //get third muon
           child = treeupsilon_part2->movePointerToTheFirstChild();
           RefCountedKinematicParticle fitMu3 = treeupsilon_part2->currentParticle();
           if(!child) break;
           TLorentzVector mu3p4;
           mu3p4.SetXYZM( fitMu3->currentState().kinematicParameters().momentum().x(), fitMu3->currentState().kinematicParameters().momentum().y(),fitMu3->currentState().kinematicParameters().momentum().z(), fitMu3->currentState().mass() );
           fourMuFit_mu3Pt.push_back(mu3p4.Pt());
           fourMuFit_mu3Eta.push_back(mu3p4.Eta());
           fourMuFit_mu3Phi.push_back(mu3p4.Phi());
           fourMuFit_mu3E.push_back(mu3p4.E());
           if (verbose) cout<<"post Fit mu3pt: "<<mu3p4.Pt()<<"mu3eta: "<<mu3p4.Eta()<<"mu3phi: "<<mu3p4.Phi()<<"mu3E: "<<mu3p4.E()<<endl;
           //get fourth muon
           child = treeupsilon_part2->movePointerToTheNextChild();
           RefCountedKinematicParticle fitMu4 = treeupsilon_part2->currentParticle();
           if(!child) break;
           TLorentzVector mu4p4;
           mu4p4.SetXYZM( fitMu4->currentState().kinematicParameters().momentum().x(), fitMu4->currentState().kinematicParameters().momentum().y(), fitMu4->currentState().kinematicParameters().momentum().z(), fitMu4->currentState().mass() );
           fourMuFit_mu4Pt.push_back(mu4p4.Pt());
           fourMuFit_mu4Eta.push_back(mu4p4.Eta());
           fourMuFit_mu4Phi.push_back(mu4p4.Phi());
           fourMuFit_mu4E.push_back(mu4p4.E());
           if (verbose) cout<<"post Fit mu4pt: "<<mu4p4.Pt()<<"mu4eta: "<<mu4p4.Eta()<<"mu4phi: "<<mu4p4.Phi()<<"mu4E: "<<mu4p4.E()<<endl;
           mu1_Pt.push_back(AllMuons[i1].pt());
           mu2_Pt.push_back(AllMuons[i2].pt());
           mu3_Pt.push_back(AllMuons[j1].pt());
           mu4_Pt.push_back(AllMuons[j2].pt());
           mu1_Eta.push_back(AllMuons[i1].eta());
           mu2_Eta.push_back(AllMuons[i2].eta());
           mu3_Eta.push_back(AllMuons[j1].eta());
           mu4_Eta.push_back(AllMuons[j2].eta());
           mu1_Phi.push_back(AllMuons[i1].phi());
           mu2_Phi.push_back(AllMuons[i2].phi());
           mu3_Phi.push_back(AllMuons[j1].phi());
           mu4_Phi.push_back(AllMuons[j2].phi());
           mu1_E.push_back(AllMuons[i1].energy());
           mu2_E.push_back(AllMuons[i2].energy());
           mu3_E.push_back(AllMuons[j1].energy());
           mu4_E.push_back(AllMuons[j2].energy());
           mu1Charge.push_back(AllMuons[i1].charge());
           mu2Charge.push_back(AllMuons[i2].charge());
           mu3Charge.push_back(AllMuons[j1].charge());
           mu4Charge.push_back(AllMuons[j2].charge());
           //mu1_genindex.push_back(mu_genindex[i1]);
           //mu2_genindex.push_back(mu_genindex[i2]);
           //mu3_genindex.push_back(mu_genindex[j1]);
           //mu4_genindex.push_back(mu_genindex[j2]);
           mu1_Soft.push_back(muon::isSoftMuon(AllMuons[i1], thePrimaryV));
           mu2_Soft.push_back(muon::isSoftMuon(AllMuons[i2], thePrimaryV));
           mu3_Soft.push_back(muon::isSoftMuon(AllMuons[j1], thePrimaryV));
           mu4_Soft.push_back(muon::isSoftMuon(AllMuons[j2], thePrimaryV));
           mu1_Tight.push_back(muon::isTightMuon(AllMuons[i1], thePrimaryV));
           mu2_Tight.push_back(muon::isTightMuon(AllMuons[i2], thePrimaryV));
           mu3_Tight.push_back(muon::isTightMuon(AllMuons[j1], thePrimaryV));
           mu4_Tight.push_back(muon::isTightMuon(AllMuons[j2], thePrimaryV));
           mu1_Medium.push_back(muon::isMediumMuon(AllMuons[i1]));
           mu2_Medium.push_back(muon::isMediumMuon(AllMuons[i2]));
           mu3_Medium.push_back(muon::isMediumMuon(AllMuons[j1]));
           mu4_Medium.push_back(muon::isMediumMuon(AllMuons[j2]));
           mu1_Loose.push_back(muon::isLooseMuon(AllMuons[i1]));
           mu2_Loose.push_back(muon::isLooseMuon(AllMuons[i2]));
           mu3_Loose.push_back(muon::isLooseMuon(AllMuons[j1]));
           mu4_Loose.push_back(muon::isLooseMuon(AllMuons[j2]));
           // Branches for check overlaping muons
           fourMuFit_mu12overlap.push_back( muon::overlap(AllMuons[i1], AllMuons[i2]) );
           fourMuFit_mu13overlap.push_back( muon::overlap(AllMuons[i1], AllMuons[j1]) );
           fourMuFit_mu14overlap.push_back( muon::overlap(AllMuons[i1], AllMuons[j2]) );
           fourMuFit_mu23overlap.push_back( muon::overlap(AllMuons[i2], AllMuons[j1]) );
           fourMuFit_mu24overlap.push_back( muon::overlap(AllMuons[i2], AllMuons[j2]) );
           fourMuFit_mu34overlap.push_back( muon::overlap(AllMuons[j1], AllMuons[j2]) );
           fourMuFit_mu12SharedSeg.push_back( muon::sharedSegments(AllMuons[i1], AllMuons[i2]) );
           fourMuFit_mu13SharedSeg.push_back( muon::sharedSegments(AllMuons[i1], AllMuons[j1]) );
           fourMuFit_mu14SharedSeg.push_back( muon::sharedSegments(AllMuons[i1], AllMuons[j2]) );
           fourMuFit_mu23SharedSeg.push_back( muon::sharedSegments(AllMuons[i2], AllMuons[j1]) );
           fourMuFit_mu24SharedSeg.push_back( muon::sharedSegments(AllMuons[i2], AllMuons[j2]) );
           fourMuFit_mu34SharedSeg.push_back( muon::sharedSegments(AllMuons[j1], AllMuons[j2])) ;
           // detach significance 
           if (verbose) cout<<"Checking PV"<<endl;
           double Vtxfrac = 0.0;
           reco::Vertex bestPV_Bang;
           Double_t PV_bestBang_X_temp = -999;
           Double_t PV_bestBang_Y_temp = -999;
           Double_t PV_bestBang_Z_temp = -999;
           Double_t PV_bestBang_XE_temp = -999;
           Double_t PV_bestBang_YE_temp = -999;
           Double_t PV_bestBang_ZE_temp = -999;
           Double_t PV_bestBang_XYE_temp= -999;
           Double_t PV_bestBang_XZE_temp= -999;
           Double_t PV_bestBang_YZE_temp= -999;
           Double_t PV_bestBang_CL_temp = -999;
           Double_t lip= 999;
           for (unsigned int i=0; i<primaryVertices_handle->size(); i++) { 
           const Vertex &PVtxBeSp = (*primaryVertices_handle)[i];
           Double_t dx = FourMuDecayVertex->position().x() - PVtxBeSp.x();
           Double_t dy = FourMuDecayVertex->position().y() - PVtxBeSp.y();
           Double_t dz = FourMuDecayVertex->position().z() - PVtxBeSp.z();
           Double_t dist = sqrt(dx*dx + dy*dy + dz*dz);
           if (verbose) cout<<"PV: X, "<<PVtxBeSp.x()<<" y: "<<PVtxBeSp.y()<<" z: "<<PVtxBeSp.z()<<endl;
           if (verbose) cout<<"PV dist: " <<dist<<endl;
           if(dist < lip)
            {
            lip = dist ;
            //Vtxfrac = Vtxfrac + 1;
            PV_bestBang_X_temp = PVtxBeSp.x();
            PV_bestBang_Y_temp = PVtxBeSp.y();
            PV_bestBang_Z_temp = PVtxBeSp.z();
            PV_bestBang_XE_temp = PVtxBeSp.covariance(0, 0);
            PV_bestBang_YE_temp = PVtxBeSp.covariance(1, 1);
            PV_bestBang_ZE_temp = PVtxBeSp.covariance(2, 2);
            PV_bestBang_XYE_temp = PVtxBeSp.covariance(0, 1);
            PV_bestBang_XZE_temp = PVtxBeSp.covariance(0, 2);
            PV_bestBang_YZE_temp = PVtxBeSp.covariance(1, 2);
            PV_bestBang_CL_temp = (TMath::Prob(PVtxBeSp.chi2(),(int)PVtxBeSp.ndof()) );
            bestPV_Bang = PVtxBeSp;
              }
             }// loop over all PV 
            vector<reco::TransientTrack> vertexTracks;
            int mu1frompv, mu2frompv, mu3frompv, mu4frompv;
            mu1frompv = 0;
            mu2frompv = 0;
            mu3frompv = 0;
            mu4frompv = 0; 
            for ( std::vector<TrackBaseRef >::const_iterator iTrack = bestPV_Bang.tracks_begin(); iTrack != bestPV_Bang.tracks_end(); ++iTrack)
              {
              TrackRef trackRef = iTrack->castTo<TrackRef>();
              TLorentzVector p4tr_RF;
              p4tr_RF.SetPtEtaPhiM(trackRef->pt(), trackRef->eta(), trackRef->phi(), muonMass);
              if ( deltaR(trackRef->eta(),trackRef->phi(),mu1p4.Eta(),mu1p4.Eta()) < 0.01 && abs(trackRef->pt() - mu1p4.Pt() ) < 0.01) mu1frompv = 1;
              if ( deltaR(trackRef->eta(),trackRef->phi(),mu2p4.Eta(),mu2p4.Eta()) < 0.01 && abs(trackRef->pt() - mu2p4.Pt() ) < 0.01) mu2frompv = 1;
              if ( deltaR(trackRef->eta(),trackRef->phi(),mu3p4.Eta(),mu3p4.Eta()) < 0.01 && abs(trackRef->pt() - mu3p4.Pt() ) < 0.01) mu3frompv = 1;
              if ( deltaR(trackRef->eta(),trackRef->phi(),mu4p4.Eta(),mu4p4.Eta()) < 0.01 && abs(trackRef->pt() - mu4p4.Pt() ) < 0.01) mu4frompv = 1;
               }  // loop over all PV 
            fourMuFit_mu1frompv.push_back(mu1frompv);
            fourMuFit_mu2frompv.push_back(mu2frompv);    
            fourMuFit_mu3frompv.push_back(mu3frompv);
            fourMuFit_mu4frompv.push_back(mu4frompv);
            fourMuFit_PVX.push_back(PV_bestBang_X_temp);
            fourMuFit_PVY.push_back(PV_bestBang_Y_temp);
            fourMuFit_PVZ.push_back(PV_bestBang_Z_temp);
            fourMuFit_PVXE.push_back(PV_bestBang_XE_temp);
            fourMuFit_PVYE.push_back(PV_bestBang_YE_temp);
            fourMuFit_PVZE.push_back(PV_bestBang_ZE_temp);
            fourMuFit_PV_XYE.push_back(PV_bestBang_XYE_temp);
            fourMuFit_PV_XZE.push_back(PV_bestBang_XZE_temp);
            fourMuFit_PV_YZE.push_back(PV_bestBang_YZE_temp);
            fourMuFit_PV_CL.push_back(PV_bestBang_CL_temp);
            if (!((trigger&16)==0))
             {
            if (verbose) cout<<"Event passed trigger selection"<<endl;
            fourmu_mass_trigger = fitFourMu->currentState().mass();
              }
           if (verbose && (muon::overlap(AllMuons[i1], AllMuons[i2])==0) && (muon::overlap(AllMuons[i1], AllMuons[j1])) && (muon::overlap(AllMuons[i1], AllMuons[j2])) && (muon::overlap(AllMuons[i2], AllMuons[j1])) && (muon::overlap(AllMuons[i2], AllMuons[j2])) && (muon::overlap(AllMuons[j1], AllMuons[j2]))) cout<<"No overlaping muons:"<<endl;
           if (verbose && !(fabs(a_wrongcomb_mass - omega_mass) < 2.0*a_wrongcomb_mass_Error) && !(fabs(b_wrongcomb_mass - omega_mass) < 2.0*b_wrongcomb_mass_Error)) cout<<"No omega in unselected pair "<<endl;
           if (verbose && !(fabs(a_wrongcomb_mass - phi_mass) < 2.0*a_wrongcomb_mass_Error) && !(fabs(b_wrongcomb_mass - phi_mass) < 2.0*b_wrongcomb_mass_Error))  cout<<"No Phi in  unselected pair "<<endl;
           if (verbose && v_mumufit_VtxCL_noMC[i] > upsilon_vtx_cut && v_mumufit_VtxCL_noMC[j] > upsilon_vtx_cut ) cout << "Passed dimuon vertex prob. cut "<<endl;
           if (verbose && temp_fourMuFit_VtxProb > 0.01) cout<<"Passed four muon vertex prob. cut "<<endl;

           if ( (muon::overlap(AllMuons[i1], AllMuons[i2])==0) && (muon::overlap(AllMuons[i1], AllMuons[j1])==0) && (muon::overlap(AllMuons[i1], AllMuons[j2])==0) && (muon::overlap(AllMuons[i2], AllMuons[j1])==0) && (muon::overlap(AllMuons[i2], AllMuons[j2])==0) && (muon::overlap(AllMuons[j1], AllMuons[j2])==0) && !(fabs(a_wrongcomb_mass - omega_mass) < 2.0*a_wrongcomb_mass_Error)&& !(fabs(b_wrongcomb_mass - omega_mass) < 2.0*b_wrongcomb_mass_Error) && !(fabs(a_wrongcomb_mass - phi_mass) < 2.0*a_wrongcomb_mass_Error) && !(fabs(b_wrongcomb_mass - phi_mass) < 2.0*b_wrongcomb_mass_Error) && temp_fourMuFit_VtxProb > 0.01 && v_mumufit_VtxCL_noMC[i] > upsilon_vtx_cut && v_mumufit_VtxCL_noMC[j] > upsilon_vtx_cut && v_mumu_Pt[i] > 3.5 && v_mumu_Pt[j] > 3.5)        
         {
         if (verbose) cout<<"Event passed all event selection"<<endl;
         fourmu_mass_allcuts = fitFourMu->currentState().mass(); 
	} // applying all cuts for four muon related reco. efficiencies 
        }// Loop2 over UPS candidates
    } // Loop1 over UPS candidates
} //end Stand Alone YY Fucntion
void MuMuGammaRootupler::fourMuonFit(pat::CompositeCandidate dimuonCand, RefCountedKinematicTree treeupsilon_part1, edm::Handle< edm::View<pat::Muon> > muons, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV){
	std::vector<pat::Muon> theRestMuons;
        vector<RefCountedKinematicParticle> Chi_1;
        treeupsilon_part1->movePointerToTheTop();
        RefCountedKinematicParticle upsilon_part1 = treeupsilon_part1->currentParticle();
        Chi_1.push_back(upsilon_part1);
        if (verbose) cout<<"Muon 1 pT: "<< dimuonCand.daughter("muon1")->pt()<< "Eta: "<<dimuonCand.daughter("muon1")->eta()<<"charge: "<<dimuonCand.daughter("muon1")->charge()<<endl;
        if (verbose) cout<<"Muon 2 pT: "<< dimuonCand.daughter("muon2")->pt()<< "Eta: "<<dimuonCand.daughter("muon2")->eta()<<"charge: "<<dimuonCand.daughter("muon2")->charge()<<endl; 
	for (edm::View<pat::Muon>::const_iterator mu3 = muons->begin(), muend = muons->end(); mu3 != muend; ++mu3){
		if (mu3->pt()<2.0 || fabs(mu3->eta())>2.4)  continue;
		if (mu3-muons->begin() == dimuonCand.userInt("mu1Index"))  continue;
		if (mu3-muons->begin() == dimuonCand.userInt("mu2Index"))  continue;
		reco::GenParticleRef genMu3;
		if (isMC_) genMu3 = mu3->genParticleRef();
		if (!isMC_ || (isMC_ && !genMu3.isNonnull())) theRestMuons.push_back(*mu3);
		for(edm::View<pat::Muon>::const_iterator mu4 = mu3+1 ; mu4 != muend; ++mu4){
			if (mu4->pt()<2.0 || fabs(mu4->eta())>2.4)  continue;
			if (mu4-muons->begin() == dimuonCand.userInt("mu1Index")) continue;
			if (mu4-muons->begin() == dimuonCand.userInt("mu2Index")) continue;
                        if (verbose) cout<<"Muon 3 pT: "<< mu3->pt()<< "Eta: "<<mu3->eta()<<endl;
                        if (verbose) cout<<"Muon 4 pT: "<< mu4->pt()<< "Eta: "<<mu4->eta()<<endl;
			reco::GenParticleRef genMu4;
			if (isMC_) genMu4 = mu4->genParticleRef();
			if (mu3->charge() == mu4->charge()) continue;
			if (!isMC_ || (isMC_ && !genMu3.isNonnull() && !genMu4.isNonnull())) {
				TLorentzVector mu1p4,mu2p4,mu3p4,mu4p4,fourMup4;
				reco::Candidate::LorentzVector v1 = dimuonCand.daughter("muon1")->p4();
				reco::Candidate::LorentzVector v2 = dimuonCand.daughter("muon2")->p4();
				mu1p4.SetPtEtaPhiM(v1.pt(),v1.eta(),v1.phi(),v1.mass());
				mu2p4.SetPtEtaPhiM(v2.pt(),v2.eta(),v2.phi(),v2.mass());
				mu3p4.SetPtEtaPhiM(mu3->pt(), mu3->eta(), mu3->phi(), mu3->mass());
				mu4p4.SetPtEtaPhiM(mu4->pt(), mu4->eta(), mu4->phi(), mu4->mass());  
                                if (isTriggerMatch_)
                                { 
                                bool Rest_Muon_trigger_Matched = TriggerMatch_restMuons(mu3p4,mu4p4);
                                if (!Rest_Muon_trigger_Matched) continue;
                                  }
				fourMup4 = mu1p4 + mu2p4 + mu3p4 + mu4p4;
				fourMuFit_Mass_allComb.push_back(fourMup4.M());
			          }

			//std::cout<<"found good mu3mu4: "<<mu3->pt()<<" "<<mu4->pt()<<", mu1: "<<muon1TT.track().pt()<<", mu2: "<<muon2TT.track().pt()<<std::endl;
			reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
			reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
			reco::TrackRef muTrack3_ref = mu3->track();
			reco::TrackRef muTrack4_ref = mu4->track();
			reco::TransientTrack muon1TT(muTrack1_ref, &(*bFieldHandle));
			reco::TransientTrack muon2TT(muTrack2_ref, &(*bFieldHandle));
			reco::TransientTrack muon3TT(muTrack3_ref, &(*bFieldHandle));
			reco::TransientTrack muon4TT(muTrack4_ref, &(*bFieldHandle));
			KinematicParticleFactoryFromTransientTrack fourMuFactory;
                        KinematicParticleFactoryFromTransientTrack smumuFactory;
			std::vector<RefCountedKinematicParticle> fourMuParticles;
                        std::vector<RefCountedKinematicParticle> smumuParticles;

                        smumuParticles.push_back(smumuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
                        smumuParticles.push_back(smumuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));


                        KinematicParticleVertexFitter smumufitter;
                        RefCountedKinematicTree smumuVertexFitTree;
                        KinematicParticleFitter csFitter;
                        smumuVertexFitTree = smumufitter.fit(smumuParticles);
                        if (!(smumuVertexFitTree->isValid())) continue;
                        if (verbose) cout<<" second dimuonCand have valid vertex"<<endl;
                        RefCountedKinematicParticle smumu_vFit_noMC;
                        RefCountedKinematicVertex smumu_vFit_vertex_noMC;
                        ParticleMass ups_m = upsilon_mass_;
                        float ups_m_sigma = upsilon_sigma_;
                        KinematicConstraint *ups_c = new MassKinematicConstraint(ups_m,ups_m_sigma);
                        smumuVertexFitTree = csFitter.fit(ups_c,smumuVertexFitTree);
                        if (!(smumuVertexFitTree->isValid())) continue;
                        if (verbose) cout<<" second dimuonCand have valid vertex after mass constrained Fit"<<endl;
                        smumuVertexFitTree->movePointerToTheTop();
                        smumu_vFit_noMC = smumuVertexFitTree->currentParticle();
                        Chi_1.push_back(smumu_vFit_noMC);
                        smumu_vFit_vertex_noMC = smumuVertexFitTree->currentDecayVertex();
                        float thisdimuon_mass = smumu_vFit_noMC->currentState().mass();
                        if (verbose) cout << "Second dimuon mass:" <<thisdimuon_mass<<endl;
                        float thisdimuon_vtxprob = ChiSquaredProbability((double)(smumu_vFit_vertex_noMC->chiSquared()),(double)(smumu_vFit_vertex_noMC->degreesOfFreedom()));
                        if (verbose) cout<<"second dimuon mass constrained vetex fit:"<<thisdimuon_vtxprob<<endl;  
                        if (thisdimuon_vtxprob < upsilon_vtx_cut) continue;
			fourMuParticles.push_back(fourMuFactory.particle (muon1TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon2TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));
		        KinematicParticleVertexFitter  Fitter;
			RefCountedKinematicTree fourMuTree = Fitter.fit(fourMuParticles);
			if(!fourMuTree->isEmpty())
			{
				fourMuTree->movePointerToTheTop();
				RefCountedKinematicParticle fitFourMu = fourMuTree->currentParticle();
				RefCountedKinematicVertex FourMuDecayVertex = fourMuTree->currentDecayVertex();
				if (fitFourMu->currentState().isValid()
						// && ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())) > fourMuFit_VtxProb)
					)
					{ //Get chib
                                                if (verbose) cout<<"Four muon mass: "<<fitFourMu->currentState().mass()<<endl;
                                                if (verbose) cout<<"Four muon vertex probabilty: "<<ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom()))<<endl;
                                                fourMuTree = Fitter.fit(Chi_1); //Two body (upsilon) Fit to get fourmu mass 
                                                if(fourMuTree->isEmpty()) continue;
                                                if(!fourMuTree->isValid()) continue; 
                                                fitFourMu = fourMuTree->currentParticle();
                                                if (verbose) cout<<"Four muon mass from two jpsi constrainted Fit: "<<fitFourMu->currentState().mass()<<endl;
						fourMuFit_Mass.push_back(fitFourMu->currentState().mass());
						fourMuFit_MassErr.push_back(sqrt(fitFourMu->currentState().kinematicParametersError().matrix()(6,6)));
						TLorentzVector p4;
						p4.SetXYZM(fitFourMu->currentState().kinematicParameters().momentum().x(),fitFourMu->currentState().kinematicParameters().momentum().y(),fitFourMu->currentState().kinematicParameters().momentum().z(),fitFourMu->currentState().mass());
						fourMuFit_Pt.push_back(p4.Pt());
						fourMuFit_Eta.push_back(p4.Eta());
						fourMuFit_Phi.push_back(p4.Phi());
						fourMuFit_VtxX.push_back(FourMuDecayVertex->position().x());
						fourMuFit_VtxY.push_back(FourMuDecayVertex->position().y());
						fourMuFit_VtxZ.push_back(FourMuDecayVertex->position().z());
						fourMuFit_VtxProb.push_back(ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())));
						fourMuFit_Chi2.push_back(FourMuDecayVertex->chiSquared());
						fourMuFit_ndof.push_back(FourMuDecayVertex->degreesOfFreedom());

						//get first muon //Childs are taken from JPsi 1 and 2
						bool child = treeupsilon_part1->movePointerToTheFirstChild();
						RefCountedKinematicParticle fitMu1 = treeupsilon_part1->currentParticle();
						if(!child) break;
						p4.SetXYZM( fitMu1->currentState().kinematicParameters().momentum().x(), fitMu1->currentState().kinematicParameters().momentum().y(), fitMu1->currentState().kinematicParameters().momentum().z(), fitMu1->currentState().mass() );
						fourMuFit_mu1Pt.push_back(p4.Pt());
						fourMuFit_mu1Eta.push_back(p4.Eta());
						fourMuFit_mu1Phi.push_back(p4.Phi());
						fourMuFit_mu1E.push_back(p4.E());

                                                if (verbose) cout<<"post Fit mu1pt: "<<p4.Pt()<<"mu1eta: "<<p4.Eta()<<"mu1phi: "<<p4.Phi()<<"mu1E: "<<p4.E()<<endl;
                                
						//get second muon
						child = treeupsilon_part1->movePointerToTheNextChild();
						RefCountedKinematicParticle fitMu2 = treeupsilon_part1->currentParticle();
						if(!child) break;
						p4.SetXYZM( fitMu2->currentState().kinematicParameters().momentum().x(), fitMu2->currentState().kinematicParameters().momentum().y(), fitMu2->currentState().kinematicParameters().momentum().z(), fitMu2->currentState().mass() );
						fourMuFit_mu2Pt.push_back(p4.Pt());
						fourMuFit_mu2Eta.push_back(p4.Eta());
						fourMuFit_mu2Phi.push_back(p4.Phi());
						fourMuFit_mu2E.push_back(p4.E());
                                                if (verbose) cout<<"post Fit mu2pt: "<<p4.Pt()<<"mu2eta: "<<p4.Eta()<<"mu2phi: "<<p4.Phi()<<"mu2E: "<<p4.E()<<endl;
						//get third muon
						child = smumuVertexFitTree->movePointerToTheFirstChild();
						RefCountedKinematicParticle fitMu3 = smumuVertexFitTree->currentParticle();
						if(!child) break;
						p4.SetXYZM( fitMu3->currentState().kinematicParameters().momentum().x(), fitMu3->currentState().kinematicParameters().momentum().y(), fitMu3->currentState().kinematicParameters().momentum().z(), fitMu3->currentState().mass() );
						fourMuFit_mu3Pt.push_back(p4.Pt());
						fourMuFit_mu3Eta.push_back(p4.Eta());
						fourMuFit_mu3Phi.push_back(p4.Phi());
						fourMuFit_mu3E.push_back(p4.E());
                                                if (verbose) cout<<"post Fit mu3pt: "<<p4.Pt()<<"mu3eta: "<<p4.Eta()<<"mu3phi: "<<p4.Phi()<<"mu3E: "<<p4.E()<<endl;
						//get fourth muon
						child = smumuVertexFitTree->movePointerToTheNextChild();
						RefCountedKinematicParticle fitMu4 = smumuVertexFitTree->currentParticle();
						if(!child) break;
						p4.SetXYZM( fitMu4->currentState().kinematicParameters().momentum().x(), fitMu4->currentState().kinematicParameters().momentum().y(), fitMu4->currentState().kinematicParameters().momentum().z(), fitMu4->currentState().mass() );
						fourMuFit_mu4Pt.push_back(p4.Pt());
						fourMuFit_mu4Eta.push_back(p4.Eta());
						fourMuFit_mu4Phi.push_back(p4.Phi());
						fourMuFit_mu4E.push_back(p4.E());
                                                if (verbose) cout<<"post Fit mu4pt: "<<p4.Pt()<<"mu4eta: "<<p4.Eta()<<"mu4phi: "<<p4.Phi()<<"mu4E: "<<p4.E()<<endl;
						mu3_Pt.push_back(mu3->pt());
						mu3_Eta.push_back(mu3->eta());
						mu3_Phi.push_back(mu3->phi());
						mu3_E.push_back(mu3->energy());
						mu4_Pt.push_back(mu4->pt());
						mu4_Eta.push_back(mu4->eta());
						mu4_Phi.push_back(mu4->phi());
						mu4_E.push_back(mu4->energy());
						reco::TrackTransientTrack muon3TTT(muTrack3_ref, &(*bFieldHandle));
						reco::TrackTransientTrack muon4TTT(muTrack4_ref, &(*bFieldHandle));
						mu3_d0.push_back(-muon3TTT.dxy(bs));
						mu3_d0err.push_back(muon3TTT.d0Error());
						mu3_dz.push_back(muon3TTT.dz());
						mu3_dzerr.push_back(muon3TTT.dzError());
						mu4_d0.push_back(-muon4TTT.dxy(bs));
						mu4_d0err.push_back(muon4TTT.d0Error());
						mu4_dz.push_back(muon4TTT.dz());
						mu4_dzerr.push_back(muon4TTT.dzError());
						mu3Charge.push_back(mu3->charge());
						mu4Charge.push_back(mu4->charge());
						mu1_Tight.push_back(tightMuon(muons->begin()+dimuonCand.userInt("mu1Index"),thePrimaryV));
						mu2_Tight.push_back(tightMuon(muons->begin()+dimuonCand.userInt("mu2Index"),thePrimaryV));
						mu3_Tight.push_back(tightMuon(mu3,thePrimaryV));
						mu4_Tight.push_back(tightMuon(mu4,thePrimaryV));
						mu1_Medium.push_back(mediumMuon(muons->begin()+dimuonCand.userInt("mu1Index")));
						mu2_Medium.push_back(mediumMuon(muons->begin()+dimuonCand.userInt("mu2Index")));
						mu3_Medium.push_back(mediumMuon(mu3));
						mu4_Medium.push_back(mediumMuon(mu4));
						mu1_Loose.push_back(looseMuon(muons->begin()+dimuonCand.userInt("mu1Index")));
						mu2_Loose.push_back(looseMuon(muons->begin()+dimuonCand.userInt("mu2Index")));
						mu3_Loose.push_back(looseMuon(mu3));
						mu4_Loose.push_back(looseMuon(mu4));

						if (isMC_) {
							reco::GenParticleRef genMu1 = (muons->begin()+dimuonCand.userInt("mu1Index"))->genParticleRef();
							reco::GenParticleRef genMu2 = (muons->begin()+dimuonCand.userInt("mu2Index"))->genParticleRef();
							if (genMu1.isNonnull() ) mu1_pdgID.push_back(genMu1->pdgId());
							else mu1_pdgID.push_back(0);
							if (genMu2.isNonnull() ) mu2_pdgID.push_back(genMu2->pdgId());
							else mu2_pdgID.push_back(0);
							//if (genMu1->motherRef()==genMu2->motherRef()) std::cout<<"genMu1->motherRef()->pdgId()="<<genMu1->motherRef()->pdgId()<<std::endl;
							//else std::cout<<"genMu1->motherRef()->pdgId()="<<genMu1->motherRef()->pdgId()<<", genMu2->motherRef()->pdgId()="<<genMu2->motherRef()->pdgId()
							// <<"genMu1->motherRef()->motherRef()->pdgId()="<<genMu1->motherRef()->motherRef()->pdgId()<<", genMu2->motherRef()->motherRef()->pdgId()="<<genMu2->motherRef()->motherRef()->pdgId()<<std::endl;

							//reco::GenParticleRef genMu3 = mu3->genParticleRef();
							//reco::GenParticleRef genMu4 = mu4->genParticleRef();
							if (genMu3.isNonnull() ){
								mu3_pdgID.push_back(genMu3->pdgId());
								/*
									if (genMu3->numberOfMothers()>0){ 
									reco::GenParticleRef mom3 = genMu3->motherRef();
									if (mom3.isNonnull()) { 
									std::cout<<""<<"mom pdgID= "<<mom3->pdgId()<<std::endl;
									if (mom3==genMu1->motherRef()) std::cout<<"same mother"<<std::endl;
									}    
									else std::cout<<"mom non"<<std::endl;
									}    
									else std::cout<<"# mom = 0"<<std::endl;
									*/
							}
							else mu3_pdgID.push_back(0);
							if (genMu4.isNonnull() ) mu4_pdgID.push_back(genMu4->pdgId());
							else mu4_pdgID.push_back(0);
						}
					}
			}
		}
	}
	muons_previousEvent.push_back(theRestMuons);
}


int MuMuGammaRootupler::fourMuonMixFit_bestYMass(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, std::vector<pat::Muon> muons_previous, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV){
	int nGoodFourMuonMix=0;
	for (std::vector<pat::Muon>::iterator mu3 = muons_previous.begin(), muend = muons_previous.end(); mu3 != muend; ++mu3){
		for(std::vector<pat::Muon>::iterator mu4 = mu3+1 ; mu4 != muend; ++mu4){
			if (mu3->charge() == mu4->charge()) continue;

			TLorentzVector mu1p4,mu2p4,mu3p4,mu4p4,fourMup4;
			reco::Candidate::LorentzVector v1 = dimuonCand.daughter("muon1")->p4();
			reco::Candidate::LorentzVector v2 = dimuonCand.daughter("muon2")->p4();
			mu1p4.SetPtEtaPhiM(v1.pt(),v1.eta(),v1.phi(),v1.mass());
			mu2p4.SetPtEtaPhiM(v2.pt(),v2.eta(),v2.phi(),v2.mass());
			mu3p4.SetPtEtaPhiM(mu3->pt(), mu3->eta(), mu3->phi(), mu3->mass());
			mu4p4.SetPtEtaPhiM(mu4->pt(), mu4->eta(), mu4->phi(), mu4->mass());
			fourMup4 = mu1p4 + mu2p4 + mu3p4 + mu4p4;
			fourMuFit_Mass_allComb_mix_bestYMass.push_back(fourMup4.M());
		}
	}
	//std::cout<<"previousSize="<<muons_previous.size()<<", thisEventSize="<<muons->size()<<std::endl;
	int muons_previous_originalSize=muons_previous.size();

	for (edm::View<pat::Muon>::const_iterator mu = muons->begin(); mu !=  muons->end(); ++mu){
		if (mu->pt()<2.0 || fabs(mu->eta())>2.4) continue;
		if (mu-muons->begin() == dimuonCand.userInt("mu1Index"))  continue;
		if (mu-muons->begin() == dimuonCand.userInt("mu2Index"))  continue;
		reco::GenParticleRef genMu;
		if (isMC_) genMu = mu->genParticleRef();
		if (!isMC_ || (isMC_ && !genMu.isNonnull())) {
			muons_previous.push_back(*mu);
		}
	}

	fourMuFit_VtxProb_mix_bestYMass = -1;
	//std::cout<<"combinedSize="<<muons_previous.size()<<std::endl;
	for (std::vector<pat::Muon>::iterator mu3 = muons_previous.begin(), muend = muons_previous.end(); mu3-muons_previous.begin() < muons_previous_originalSize; ++mu3){
		for(std::vector<pat::Muon>::iterator mu4 = mu3+1 ; mu4 != muend; ++mu4){
			if (mu3->charge() == mu4->charge()) continue;

			reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
			reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
			reco::TrackRef muTrack3_ref = mu3->track();
			reco::TrackRef muTrack4_ref = mu4->track();
			reco::TransientTrack muon1TT(muTrack1_ref, &(*bFieldHandle));
			reco::TransientTrack muon2TT(muTrack2_ref, &(*bFieldHandle));
			reco::TransientTrack muon3TT(muTrack3_ref, &(*bFieldHandle));
			reco::TransientTrack muon4TT(muTrack4_ref, &(*bFieldHandle));
			KinematicParticleFactoryFromTransientTrack fourMuFactory;
			std::vector<RefCountedKinematicParticle> fourMuParticles;
			fourMuParticles.push_back(fourMuFactory.particle (muon1TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon2TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));
			KinematicConstrainedVertexFitter constVertexFitter;
			//fit w/ mass constraint
			//MultiTrackKinematicConstraint *upsilon_mtc = new  TwoTrackMassKinematicConstraint(upsilon_mass_);
			//RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles,upsilon_mtc);

			//fit w/o any mass constraint
			RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles);

			if(!fourMuTree->isEmpty())
			{
				fourMuTree->movePointerToTheTop();
				RefCountedKinematicParticle fitFourMu = fourMuTree->currentParticle();
				RefCountedKinematicVertex FourMuDecayVertex = fourMuTree->currentDecayVertex();
				if (fitFourMu->currentState().isValid() &&
						ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())) > fourMuFit_VtxProb_mix_bestYMass)
				{ //Get chib         
					nGoodFourMuonMix = 1;
					fourMuFit_Mass_mix_bestYMass = fitFourMu->currentState().mass();
					fourMuFit_MassErr_mix_bestYMass = sqrt(fitFourMu->currentState().kinematicParametersError().matrix()(6,6));
					fourMuFit_p4_mix_bestYMass.SetXYZM(fitFourMu->currentState().kinematicParameters().momentum().x(),fitFourMu->currentState().kinematicParameters().momentum().y(),fitFourMu->currentState().kinematicParameters().momentum().z(),fourMuFit_Mass_mix_bestYMass);
					fourMuFit_VtxX_mix_bestYMass = FourMuDecayVertex->position().x();
					fourMuFit_VtxY_mix_bestYMass = FourMuDecayVertex->position().y();
					fourMuFit_VtxZ_mix_bestYMass = FourMuDecayVertex->position().z();
					fourMuFit_VtxProb_mix_bestYMass = ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom()));
					fourMuFit_Chi2_mix_bestYMass = FourMuDecayVertex->chiSquared();
					fourMuFit_ndof_mix_bestYMass = FourMuDecayVertex->degreesOfFreedom();
					if (mu4-muons_previous.begin() >= muons_previous_originalSize) fourMuFit_3plus1_mix_bestYMass = 1;
					else fourMuFit_3plus1_mix_bestYMass = 0;

					//get first muon
					bool child = fourMuTree->movePointerToTheFirstChild();
					RefCountedKinematicParticle fitMu1 = fourMuTree->currentParticle();
					if(!child) break;
					double mu1M_fit_mix = fitMu1->currentState().mass();
					double mu1Px_fit_mix = fitMu1->currentState().kinematicParameters().momentum().x();
					double mu1Py_fit_mix = fitMu1->currentState().kinematicParameters().momentum().y();
					double mu1Pz_fit_mix = fitMu1->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu1p4_mix_bestYMass.SetXYZM( mu1Px_fit_mix, mu1Py_fit_mix, mu1Pz_fit_mix, mu1M_fit_mix );
					//get second muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu2 = fourMuTree->currentParticle();
					if(!child) break;
					double mu2M_fit_mix = fitMu2->currentState().mass();
					double mu2Px_fit_mix = fitMu2->currentState().kinematicParameters().momentum().x();
					double mu2Py_fit_mix = fitMu2->currentState().kinematicParameters().momentum().y();
					double mu2Pz_fit_mix = fitMu2->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu2p4_mix_bestYMass.SetXYZM( mu2Px_fit_mix, mu2Py_fit_mix, mu2Pz_fit_mix, mu2M_fit_mix );

					//get third muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu3 = fourMuTree->currentParticle();
					if(!child) break;
					double mu3M_fit_mix = fitMu3->currentState().mass();
					double mu3Px_fit_mix = fitMu3->currentState().kinematicParameters().momentum().x();
					double mu3Py_fit_mix = fitMu3->currentState().kinematicParameters().momentum().y();
					double mu3Pz_fit_mix = fitMu3->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu3p4_mix_bestYMass.SetXYZM( mu3Px_fit_mix, mu3Py_fit_mix, mu3Pz_fit_mix, mu3M_fit_mix );

					//get fourth muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu4 = fourMuTree->currentParticle();
					if(!child) break;
					double mu4M_fit_mix = fitMu4->currentState().mass();
					double mu4Px_fit_mix = fitMu4->currentState().kinematicParameters().momentum().x();
					double mu4Py_fit_mix = fitMu4->currentState().kinematicParameters().momentum().y();
					double mu4Pz_fit_mix = fitMu4->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu4p4_mix_bestYMass.SetXYZM( mu4Px_fit_mix, mu4Py_fit_mix, mu4Pz_fit_mix, mu4M_fit_mix );

					reco::Candidate::LorentzVector v3 = mu3->p4();
					reco::Candidate::LorentzVector v4 = mu4->p4();
					mu3_p4_mix_bestYMass.SetPtEtaPhiM(v3.pt(),v3.eta(),v3.phi(),v3.mass());
					mu4_p4_mix_bestYMass.SetPtEtaPhiM(v4.pt(),v4.eta(),v4.phi(),v4.mass());
					reco::TrackTransientTrack muon3TTT(muTrack3_ref, &(*bFieldHandle));
					reco::TrackTransientTrack muon4TTT(muTrack4_ref, &(*bFieldHandle));
					mu3_d0_mix_bestYMass = -muon3TTT.dxy(bs);
					mu3_d0err_mix_bestYMass = muon3TTT.d0Error();
					mu3_dz_mix_bestYMass = muon3TTT.dz();
					mu3_dzerr_mix_bestYMass = muon3TTT.dzError();
					mu4_d0_mix_bestYMass = -muon4TTT.dxy(bs);
					mu4_d0err_mix_bestYMass = muon4TTT.d0Error();
					mu4_dz_mix_bestYMass = muon4TTT.dz();
					mu4_dzerr_mix_bestYMass = muon4TTT.dzError();
					mu3Charge_mix_bestYMass = mu3->charge();
					mu4Charge_mix_bestYMass = mu4->charge();
				}
			}
		}
	}
	return nGoodFourMuonMix;
}


void MuMuGammaRootupler::fourMuonFit_bestYMass(pat::CompositeCandidate dimuonCand, RefCountedKinematicTree treeupsilon_part1, edm::Handle< edm::View<pat::Muon> > muons, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV){
	muons_previousEvent_bestYMass.clear();
	fourMuFit_VtxProb_bestYMass = -1;
        vector<RefCountedKinematicParticle> Chi_1;
        treeupsilon_part1->movePointerToTheTop();
        RefCountedKinematicParticle upsilon_part1 = treeupsilon_part1->currentParticle();
        Chi_1.push_back(upsilon_part1);

	//std::cout<<"mu1Index="<<mu1Index<<", mu2Index="<<mu2Index<<std::endl; 
	for (edm::View<pat::Muon>::const_iterator mu3 = muons->begin(), muend = muons->end(); mu3 != muend; ++mu3){
		if (mu3->pt()<2.0 || fabs(mu3->eta())>2.4)  continue;
		if (mu3-muons->begin() == dimuonCand.userInt("mu1Index"))  continue;
		if (mu3-muons->begin() == dimuonCand.userInt("mu2Index"))  continue;
		reco::GenParticleRef genMu3;
		if (isMC_) genMu3 = mu3->genParticleRef();
		if (!isMC_ || (isMC_ && !genMu3.isNonnull())) muons_previousEvent_bestYMass.push_back(*mu3);
		//std::cout<<"fill vector: "<<mu3->pt()<<std::endl;
		for(edm::View<pat::Muon>::const_iterator mu4 = mu3+1 ; mu4 != muend; ++mu4){
			if (mu4->pt()<2.0 || fabs(mu4->eta())>2.4)  continue;
			if (mu4-muons->begin() == dimuonCand.userInt("mu1Index")) continue;
			if (mu4-muons->begin() == dimuonCand.userInt("mu2Index")) continue;
			reco::GenParticleRef genMu4;
			if (isMC_) genMu4 = mu4->genParticleRef();

			if (mu3->charge() == mu4->charge()) continue;
			/*if ( (tightMuon(muons->begin()+mu1Index)+
			  tightMuon(muons->begin()+mu2Index)+
			  tightMuon(mu3)+
			  tightMuon(mu4)) < 2 
			  ) continue;*/

			if (!isMC_ || (isMC_ && !genMu3.isNonnull() && !genMu4.isNonnull())) {
				TLorentzVector mu1p4,mu2p4,mu3p4,mu4p4,fourMup4;
				reco::Candidate::LorentzVector v1 = dimuonCand.daughter("muon1")->p4();
				reco::Candidate::LorentzVector v2 = dimuonCand.daughter("muon2")->p4();
				mu1p4.SetPtEtaPhiM(v1.pt(),v1.eta(),v1.phi(),v1.mass());
				mu2p4.SetPtEtaPhiM(v2.pt(),v2.eta(),v2.phi(),v2.mass());
				mu3p4.SetPtEtaPhiM(mu3->pt(), mu3->eta(), mu3->phi(), mu3->mass());
				mu4p4.SetPtEtaPhiM(mu4->pt(), mu4->eta(), mu4->phi(), mu4->mass());
                                if (isTriggerMatch_)
                                {
                                bool Rest_Muon_trigger_Matched = TriggerMatch_restMuons(mu3p4,mu4p4);
                                if (!Rest_Muon_trigger_Matched ) continue;
                                   }
				fourMup4 = mu1p4 + mu2p4 + mu3p4 + mu4p4;
				fourMuFit_Mass_allComb_bestYMass.push_back(fourMup4.M());
			}

			//std::cout<<"found good mu3mu4: "<<mu3->pt()<<" "<<mu4->pt()<<", mu1: "<<muon1TT.track().pt()<<", mu2: "<<muon2TT.track().pt()<<std::endl;
			reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
			reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
			reco::TrackRef muTrack3_ref = mu3->track();
			reco::TrackRef muTrack4_ref = mu4->track();
			reco::TransientTrack muon1TT(muTrack1_ref, &(*bFieldHandle));
			reco::TransientTrack muon2TT(muTrack2_ref, &(*bFieldHandle));
			reco::TransientTrack muon3TT(muTrack3_ref, &(*bFieldHandle));
			reco::TransientTrack muon4TT(muTrack4_ref, &(*bFieldHandle));
			KinematicParticleFactoryFromTransientTrack fourMuFactory;
                        KinematicParticleFactoryFromTransientTrack smumuFactory;
			std::vector<RefCountedKinematicParticle> fourMuParticles;
                        std::vector<RefCountedKinematicParticle> smumuParticles;
       
                        smumuParticles.push_back(smumuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
                        smumuParticles.push_back(smumuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));

                        KinematicParticleVertexFitter smumufitter;
                        RefCountedKinematicTree smumuVertexFitTree;
                        KinematicParticleFitter csFitter;
                        smumuVertexFitTree = smumufitter.fit(smumuParticles);
                        if (!(smumuVertexFitTree->isValid())) continue;
                        if (verbose) cout<<" second dimuonCand have valid vertex"<<endl;
                        RefCountedKinematicParticle smumu_vFit_noMC;
                        RefCountedKinematicVertex smumu_vFit_vertex_noMC;
                        ParticleMass ups_m = upsilon_mass_;
                        float ups_m_sigma = upsilon_sigma_;
                        KinematicConstraint *ups_c = new MassKinematicConstraint(ups_m,ups_m_sigma);
                        smumuVertexFitTree = csFitter.fit(ups_c,smumuVertexFitTree);
                        smumuVertexFitTree->movePointerToTheTop();
                        smumu_vFit_noMC = smumuVertexFitTree->currentParticle();
                        Chi_1.push_back(smumu_vFit_noMC);
                        smumu_vFit_vertex_noMC = smumuVertexFitTree->currentDecayVertex();
                        float thisdimuon_mass = smumu_vFit_noMC->currentState().mass();
                        if (verbose) cout << "Second dimuon mass:" <<thisdimuon_mass<<endl;
                        float thisdimuon_vtxprob = ChiSquaredProbability((double)(smumu_vFit_vertex_noMC->chiSquared()),(double)(smumu_vFit_vertex_noMC->degreesOfFreedom()));
                        if (verbose) cout<<"second dimuon vetex fit:"<<thisdimuon_vtxprob<<endl;
                        if (thisdimuon_vtxprob < upsilon_vtx_cut) continue;
                        
			fourMuParticles.push_back(fourMuFactory.particle (muon1TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon2TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));

                        KinematicParticleVertexFitter  Fitter;
                        RefCountedKinematicTree fourMuTree = Fitter.fit(fourMuParticles);

			if(!fourMuTree->isEmpty())
			{
				fourMuTree->movePointerToTheTop();
				RefCountedKinematicParticle fitFourMu = fourMuTree->currentParticle();
				RefCountedKinematicVertex FourMuDecayVertex = fourMuTree->currentDecayVertex();

				if (fitFourMu->currentState().isValid())
						// && ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())) > fourMuFit_VtxProb_bestYMass)
				{ //Get chib
                                                if (verbose) cout<<"Four muon mass: "<<fitFourMu->currentState().mass()<<endl;
                                                if (verbose) cout<<"Four muon vertex probabilty"<<ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom()))<<endl;
                                                fourMuTree = Fitter.fit(Chi_1);
                                                if(fourMuTree->isEmpty()) continue;
                                                if(!fourMuTree->isValid()) continue;
                                                fitFourMu = fourMuTree->currentParticle();
                                                if (verbose) cout<<"Four muon mass from two jpsi constrainted Fit: "<<fitFourMu->currentState().mass()<<endl;

					fourMuFit_Mass_bestYMass = fitFourMu->currentState().mass();
					fourMuFit_MassErr_bestYMass = sqrt(fitFourMu->currentState().kinematicParametersError().matrix()(6,6));
					fourMuFit_p4_bestYMass.SetXYZM(fitFourMu->currentState().kinematicParameters().momentum().x(),fitFourMu->currentState().kinematicParameters().momentum().y(),fitFourMu->currentState().kinematicParameters().momentum().z(),fourMuFit_Mass_bestYMass);
					fourMuFit_VtxX_bestYMass = FourMuDecayVertex->position().x();
					fourMuFit_VtxY_bestYMass = FourMuDecayVertex->position().y();
					fourMuFit_VtxZ_bestYMass = FourMuDecayVertex->position().z();
					fourMuFit_VtxProb_bestYMass = ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom()));
					fourMuFit_Chi2_bestYMass = FourMuDecayVertex->chiSquared();
					fourMuFit_ndof_bestYMass = FourMuDecayVertex->degreesOfFreedom();

					//get first muon
					bool child = treeupsilon_part1->movePointerToTheFirstChild();
					RefCountedKinematicParticle fitMu1 = treeupsilon_part1->currentParticle();
					if(!child) break;
					double mu1M_fit = fitMu1->currentState().mass();
					double mu1Px_fit = fitMu1->currentState().kinematicParameters().momentum().x();
					double mu1Py_fit = fitMu1->currentState().kinematicParameters().momentum().y();
					double mu1Pz_fit = fitMu1->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu1p4_bestYMass.SetXYZM( mu1Px_fit, mu1Py_fit, mu1Pz_fit, mu1M_fit );

					//get second muon
					child = treeupsilon_part1->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu2 = treeupsilon_part1->currentParticle();
					if(!child) break;
					double mu2M_fit = fitMu2->currentState().mass();
					double mu2Px_fit = fitMu2->currentState().kinematicParameters().momentum().x();
					double mu2Py_fit = fitMu2->currentState().kinematicParameters().momentum().y();
					double mu2Pz_fit = fitMu2->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu2p4_bestYMass.SetXYZM( mu2Px_fit, mu2Py_fit, mu2Pz_fit, mu2M_fit );

					//get third muon
					child = smumuVertexFitTree->movePointerToTheFirstChild();
					RefCountedKinematicParticle fitMu3 = smumuVertexFitTree->currentParticle();
					if(!child) break;
					double mu3M_fit = fitMu3->currentState().mass();
					double mu3Px_fit = fitMu3->currentState().kinematicParameters().momentum().x();
					double mu3Py_fit = fitMu3->currentState().kinematicParameters().momentum().y();
					double mu3Pz_fit = fitMu3->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu3p4_bestYMass.SetXYZM( mu3Px_fit, mu3Py_fit, mu3Pz_fit, mu3M_fit );

					//get fourth muon
					child = smumuVertexFitTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu4 = smumuVertexFitTree->currentParticle();
					if(!child) break;
					double mu4M_fit = fitMu4->currentState().mass();
					double mu4Px_fit = fitMu4->currentState().kinematicParameters().momentum().x();
					double mu4Py_fit = fitMu4->currentState().kinematicParameters().momentum().y();
					double mu4Pz_fit = fitMu4->currentState().kinematicParameters().momentum().z();
					fourMuFit_mu4p4_bestYMass.SetXYZM( mu4Px_fit, mu4Py_fit, mu4Pz_fit, mu4M_fit );

					//std::cout<<fourMuFit_Mass<<" "<<(fourMuFit_mu1p4+fourMuFit_mu2p4+fourMuFit_mu3p4+fourMuFit_mu4p4).M()<<std::endl;

					reco::Candidate::LorentzVector v3 = mu3->p4();
					reco::Candidate::LorentzVector v4 = mu4->p4();
					mu3_p4_bestYMass.SetPtEtaPhiM(v3.pt(),v3.eta(),v3.phi(),v3.mass());
					mu4_p4_bestYMass.SetPtEtaPhiM(v4.pt(),v4.eta(),v4.phi(),v4.mass());
					reco::TrackTransientTrack muon3TTT(muTrack3_ref, &(*bFieldHandle));
					reco::TrackTransientTrack muon4TTT(muTrack4_ref, &(*bFieldHandle));
					mu3_d0_bestYMass = -muon3TTT.dxy(bs);
					mu3_d0err_bestYMass = muon3TTT.d0Error();
					mu3_dz_bestYMass = muon3TTT.dz();
					mu3_dzerr_bestYMass = muon3TTT.dzError();
					mu4_d0_bestYMass = -muon4TTT.dxy(bs);
					mu4_d0err_bestYMass = muon4TTT.d0Error();
					mu4_dz_bestYMass = muon4TTT.dz();
					mu4_dzerr_bestYMass = muon4TTT.dzError();
					mu3Charge_bestYMass = mu3->charge();
					mu4Charge_bestYMass = mu4->charge();
					mu1_Tight_bestYMass = tightMuon(muons->begin()+dimuonCand.userInt("mu1Index"),thePrimaryV);
					mu2_Tight_bestYMass = tightMuon(muons->begin()+dimuonCand.userInt("mu2Index"),thePrimaryV);
					mu3_Tight_bestYMass = tightMuon(mu3,thePrimaryV);
					mu4_Tight_bestYMass = tightMuon(mu4,thePrimaryV);
                                        mu1_Medium_bestYMass = mediumMuon(muons->begin()+dimuonCand.userInt("mu1Index"));
                                        mu2_Medium_bestYMass = mediumMuon(muons->begin()+dimuonCand.userInt("mu2Index"));
                                        mu3_Medium_bestYMass = mediumMuon(mu3);
                                        mu4_Medium_bestYMass = mediumMuon(mu4);
                                        mu1_Loose_bestYMass = looseMuon(muons->begin()+dimuonCand.userInt("mu1Index"));
                                        mu2_Loose_bestYMass = looseMuon(muons->begin()+dimuonCand.userInt("mu2Index"));
                                        mu3_Loose_bestYMass = looseMuon(mu3);
                                        mu4_Loose_bestYMass = looseMuon(mu4);                                        

					if (isMC_) {
						//reco::GenParticleRef genMu1 = (muons->begin()+mu1Index)->genParticleRef();
						//reco::GenParticleRef genMu2 = (muons->begin()+mu2Index)->genParticleRef();
						//if (genMu1->motherRef()==genMu2->motherRef()) std::cout<<"genMu1->motherRef()->pdgId()="<<genMu1->motherRef()->pdgId()<<std::endl;
						//else std::cout<<"genMu1->motherRef()->pdgId()="<<genMu1->motherRef()->pdgId()<<", genMu2->motherRef()->pdgId()="<<genMu2->motherRef()->pdgId()
						// <<"genMu1->motherRef()->motherRef()->pdgId()="<<genMu1->motherRef()->motherRef()->pdgId()<<", genMu2->motherRef()->motherRef()->pdgId()="<<genMu2->motherRef()->motherRef()->pdgId()<<std::endl;

						//reco::GenParticleRef genMu3 = mu3->genParticleRef();
						//reco::GenParticleRef genMu4 = mu4->genParticleRef();
						if (genMu3.isNonnull() ){
							mu3_pdgID_bestYMass = genMu3->pdgId();
							        /*
								if (genMu3->numberOfMothers()>0){ 
								reco::GenParticleRef mom3 = genMu3->motherRef();
								if (mom3.isNonnull()) { 
								std::cout<<""<<"mom pdgID= "<<mom3->pdgId()<<std::endl;
								if (mom3==genMu1->motherRef()) std::cout<<"same mother"<<std::endl;
								}    
								else std::cout<<"mom non"<<std::endl;
								}    
								else std::cout<<"# mom = 0"<<std::endl;
								*/
						}
						else mu3_pdgID_bestYMass = 0;
						if (genMu4.isNonnull() ) mu4_pdgID_bestYMass = genMu4->pdgId();
						else mu4_pdgID_bestYMass = 0;
					}
				}
			}
		}
	}
}
double MuMuGammaRootupler::GetcTau(RefCountedKinematicVertex& decayVrtx, RefCountedKinematicParticle& kinePart, Vertex& bs, bool MassConstraint)
{
    TVector3 vtx;
    TVector3 pvtx;
    vtx.SetXYZ((*decayVrtx).position().x(), (*decayVrtx).position().y(), 0);
    pvtx.SetXYZ(bs.position().x(), bs.position().y(), 0);
    TVector3 pperp(kinePart->currentState().globalMomentum().x(),
        kinePart->currentState().globalMomentum().y(), 0);
    TVector3 vdiff = vtx - pvtx;

    //Jinfeng 10.3
    GlobalError DecayErr = (*decayVrtx).error();
    GlobalError PrimaryErr = bs.error();

    AlgebraicSymMatrix33 LxyErrMatrix = DecayErr.matrix() + PrimaryErr.matrix();
    AlgebraicVector3 vdiff_;
    vdiff_[0] = vdiff.x();
    vdiff_[1] = vdiff.y();
    vdiff_[2] = 0;

    double LxyErr2 = ROOT::Math::Similarity(LxyErrMatrix, vdiff_);
    double LxyPVErr = 0;
    double LxyPVSig = 0;

    double cosAlpha = vdiff.Dot(pperp) / (vdiff.Perp() * pperp.Perp());
    double LxyPV = vdiff.Dot(pperp) / pperp.Mag();
    double Lxy = vdiff.Perp();
    double ctauPV = LxyPV * kinePart->currentState().mass() / pperp.Perp();

    if (Lxy != 0)
        LxyPVErr = TMath::Sqrt(LxyErr2) / Lxy * cosAlpha;

    if (LxyPVErr != 0)
        LxyPVSig = LxyPV / LxyPVErr;

    //Jinfeng
    if (MassConstraint) {
        temp_Lxy_MC = Lxy;
        temp_LxyPV_MC = LxyPV;
        temp_CosAlpha_MC = cosAlpha;
        temp_LxyPVErr_MC = LxyPVErr;
        temp_LxyPVSig_MC = LxyPVSig;
    }
    else {
        temp_Lxy_noMC = Lxy;
        temp_LxyPV_noMC = LxyPV;
        temp_CosAlpha_noMC = cosAlpha;
        temp_LxyPVErr_noMC = LxyPVErr;
        temp_LxyPVSig_noMC = LxyPVSig;
    }

    return ctauPV;
    //return distXY.value();
  }
double MuMuGammaRootupler::GetcTauErr(RefCountedKinematicVertex& decayVrtx, RefCountedKinematicParticle& kinePart, Vertex& bs)
  {
    TVector3 pperp(kinePart->currentState().globalMomentum().x(),
                   kinePart->currentState().globalMomentum().y(), 0);
    AlgebraicVector vpperp(3);
    vpperp[0] = pperp.x();
    vpperp[1] = pperp.y();
    vpperp[2] = 0.;

    GlobalError v1e = (Vertex(*decayVrtx)).error();
    GlobalError v2e = bs.error();
    AlgebraicSymMatrix vXYe = asHepMatrix(v1e.matrix()) + asHepMatrix(v2e.matrix());
    double ctauErrPV = sqrt(vXYe.similarity(vpperp)) * kinePart->currentState().mass() / (pperp.Perp2());

    return ctauErrPV;
  }


//Jinfeng 10.3
void MuMuGammaRootupler::GetDistance(RefCountedKinematicVertex decayVrtx1, RefCountedKinematicVertex decayVrtx2, bool MassConstraint)
{
    TVector3 DistanceVector;
    DistanceVector.SetXYZ((*decayVrtx1).position().x() - (*decayVrtx2).position().x(), (*decayVrtx1).position().y() - (*decayVrtx2).position().y(), 0);

    double Distance_ = DistanceVector.Perp();

    GlobalError DecayErr1 = (*decayVrtx1).error();
    GlobalError DecayErr2 = (*decayVrtx2).error();

    AlgebraicSymMatrix33 DistanceErrMatrix = DecayErr1.matrix() + DecayErr2.matrix();
    AlgebraicVector3 DistanceVector_;
    DistanceVector_[0] = DistanceVector.x();
    DistanceVector_[1] = DistanceVector.y();
    DistanceVector_[2] = 0;

    double DistanceErr2 = ROOT::Math::Similarity(DistanceErrMatrix, DistanceVector_);
    double DistanceErr = 0;
    double DistanceSig = 0;

    if (Distance_ != 0)
        DistanceErr = TMath::Sqrt(DistanceErr2) / Distance_;

    if (DistanceErr != 0)
        DistanceSig = Distance_ / DistanceErr;

    if (MassConstraint) {
        temp_Distance_MC = Distance_;
        temp_DistanceErr_MC = DistanceErr;
        temp_DistanceSig_MC = DistanceSig;
    }

    else {
        temp_Distance_noMC = Distance_;
        temp_DistanceErr_noMC = DistanceErr;
        temp_DistanceSig_noMC = DistanceSig;
    }

    return;
}


//define this as a plug-in
DEFINE_FWK_MODULE(MuMuGammaRootupler);
