#include "ReactionProcess.hh"
#include "DetectorConstruction.hh"

#include "G4DynamicParticle.hh"
#include "G4EventManager.hh"
#include "G4IonTable.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleTable.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "Analysis.hh"

#include "TLorentzVector.h"
#include "TVector3.h"

#include <algorithm>
#include <cmath>
#include <cfloat>
#include <stdexcept>

namespace
{
    G4ParticleDefinition *GetReactionParticleDefinition(G4int z, G4int a,
                                                        G4double excitationEnergy = 0.)
    {
        auto *particleTable = G4ParticleTable::GetParticleTable();
        if (excitationEnergy == 0.)
        {
            if (z == 0 && a == 1)
                return particleTable->FindParticle("neutron");
            if (z == 1 && a == 1)
                return particleTable->FindParticle("proton");
            if (z == 1 && a == 2)
                return particleTable->FindParticle("deuteron");
            if (z == 1 && a == 3)
                return particleTable->FindParticle("triton");
            if (z == 2 && a == 3)
                return particleTable->FindParticle("He3");
            if (z == 2 && a == 4)
                return particleTable->FindParticle("alpha");
        }

        if (z == 0)
            return particleTable->FindParticle("neutron");
        return G4IonTable::GetIonTable()->GetIon(z, a, excitationEnergy);
    }

    G4int RoundedParam(const std::map<std::string, G4double> &params,
                       const std::string &name)
    {
        const auto it = params.find(name);
        G4cout << name << G4endl;
        if (it == params.end())
           throw std::runtime_error("Missing reaction parameter: " + name);
        return static_cast<G4int>(std::lround(it->second));
    }

    G4double OptionalParam(const std::map<std::string, G4double> &params,
                           const std::string &name,
                           G4double defaultValue)
    {
        const auto it = params.find(name);
        return it == params.end() ? defaultValue : it->second;
    }
}

namespace MdmSim
{
    ReactionProcess::ReactionProcess(const G4String &processName)
        : G4VDiscreteProcess(processName, fHadronic)
    {
        SetProcessSubType(111);
    }

    ReactionProcess::~ReactionProcess() = default;

    G4double ReactionProcess::GetMeanFreePath(const G4Track &aTrack,
                                              G4double,
                                              G4ForceCondition *condition)
    {
        *condition = NotForced;

        const auto *definition = aTrack.GetParticleDefinition();
        const G4bool isBeamIon =
            definition &&
            definition->GetAtomicNumber() == fBeamZ &&
            definition->GetAtomicMass() == fBeamA &&
            aTrack.GetParentID() == 0;
        const auto *volume = aTrack.GetVolume();
        const G4bool inTarget =
            volume && volume->GetLogicalVolume()->GetName() == "Target";

        if (!isBeamIon || !inTarget || !fWillReact)
            return DBL_MAX;

        const G4double pathLength = TargetPathLength(aTrack);
        if (pathLength <= 0.)
            return DBL_MAX;

        const G4double reactionPath = fReactionFraction * pathLength;
        const G4double depth = TargetDepthAlongTrack(aTrack);
        return std::max(0., reactionPath - depth);
    }

    G4VParticleChange *ReactionProcess::PostStepDoIt(const G4Track &aTrack,
                                                     const G4Step &aStep)
    {
        aParticleChange.Initialize(aTrack);
        const G4ThreeVector reactionPosition = aStep.GetPostStepPoint()->GetPosition();
        TwoBody(aTrack, reactionPosition);
        aParticleChange.ProposeTrackStatus(fStopAndKill);
        aParticleChange.ProposeEnergy(0.);
        return &aParticleChange;
    }

    void ReactionProcess::StartTracking(G4Track *track)
    {
        G4VProcess::StartTracking(track);
        const auto *definition = track ? track->GetParticleDefinition() : nullptr;
        const G4bool isBeamIon =
            definition &&
            definition->GetAtomicNumber() == fBeamZ &&
            definition->GetAtomicMass() == fBeamA &&
            track->GetParentID() == 0;
        fWillReact = isBeamIon && G4UniformRand() < fReactionProbability;
        fReactionFraction = G4UniformRand();
        // fReactionFraction = 1; //Just to get the average energy of the detected particles
    }

    void ReactionProcess::SetParams(std::map<std::string, G4double> params)
    {
        fBeamZ = RoundedParam(params, "BeamZ");
        fBeamA = RoundedParam(params, "BeamA");
        fTargetZ = RoundedParam(params, "TargetZ");
        fTargetA = RoundedParam(params, "TargetA");
        fLightProductZ = RoundedParam(params, "LightProductZ");
        fLightProductA = RoundedParam(params, "LightProductA");
        fLightProductCharge = RoundedParam(params, "LightProductCharge");
        fLightProductEx = OptionalParam(params, "LightProductExMeV", 0.) * MeV;
        fHeavyProductZ = RoundedParam(params, "HeavyProductZ");
        fHeavyProductA = RoundedParam(params, "HeavyProductA");
        fHeavyProductCharge = RoundedParam(params, "HeavyProductCharge");
        fHeavyProductEx = RoundedParam(params, "HeavyProductExMeV") * MeV;

        fDecayHeavyZ = RoundedParam(params, "DecayHeavyZ");
        fDecayQ = RoundedParam(params, "DecayQ") * MeV;
        fDecayHeavyA = RoundedParam(params, "DecayHeavyA");
        fDecayHeavyCharge = RoundedParam(params, "DecayHeavyCharge");
        fDecayLightOneZ = RoundedParam(params, "DecayLightOneZ");
        fDecayLightOneA = RoundedParam(params, "DecayLightOneA");
        fDecayLightOneCharge = RoundedParam(params, "DecayLightOneCharge");
        fDecayLightTwoZ = RoundedParam(params, "DecayLightTwoZ");
        fDecayLightTwoA = RoundedParam(params, "DecayLightTwoA");
        fDecayLightTwoCharge = RoundedParam(params, "DecayLightTwoCharge");



        fReactionProbability = OptionalParam(params, "ReactionProbability", 1.);

        if (fBeamZ <= 0 || fBeamA <= 0 || fTargetZ < 0 || fTargetA <= 0 ||
            fLightProductZ < 0 || fLightProductA <= 0 ||
            fHeavyProductZ < 0 || fHeavyProductA <= 0)
        {
            throw std::runtime_error("Reaction Z/A values are invalid");
        }
        if (fLightProductCharge < 0 || fLightProductCharge > fLightProductZ ||
            fHeavyProductCharge < 0 || fHeavyProductCharge > fHeavyProductZ)
        {
            throw std::runtime_error("Reaction product charge states must satisfy 0 <= charge <= Z");
        }
        if (fReactionProbability < 0. || fReactionProbability > 1.)
        {
            throw std::runtime_error("ReactionProbability must satisfy 0 <= ReactionProbability <= 1");
        }
    }

    G4double ReactionProcess::TargetPathLength(const G4Track &aTrack) const
    {
        const auto *detectorConstruction =
            static_cast<const DetectorConstruction *>(
                G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        const G4double targetDz = detectorConstruction->GetTargetThicknessLength();
        G4ThreeVector localDirection = aTrack.GetMomentumDirection();
        localDirection.rotateY(-detectorConstruction->GetTargetRotationAngle());
        const G4double localCosZ = std::abs(localDirection.z());
        if (localCosZ < 1e-9)
            return DBL_MAX;
        return targetDz / localCosZ;
    }

    G4double ReactionProcess::TargetDepthAlongTrack(const G4Track &aTrack) const
    {
        const auto *detectorConstruction =
            static_cast<const DetectorConstruction *>(
                G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        const G4double targetDz = detectorConstruction->GetTargetThicknessLength();

        G4ThreeVector localPosition = aTrack.GetPosition();
        G4ThreeVector localDirection = aTrack.GetMomentumDirection();
        localPosition.rotateY(-detectorConstruction->GetTargetRotationAngle());
        localDirection.rotateY(-detectorConstruction->GetTargetRotationAngle());

        const G4double localCosZ = localDirection.z();
        if (std::abs(localCosZ) < 1e-9)
            return 0.;

        const G4double entranceZ = localCosZ >= 0. ? -0.5 * targetDz : 0.5 * targetDz;
        const G4double depth = (localPosition.z() - entranceZ) / localCosZ;
        return std::clamp(depth, 0., TargetPathLength(aTrack));
    }

    void ReactionProcess::TwoBody(const G4Track &aTrack, const G4ThreeVector &reactionPosition)
    {
        Analysis *analysis = Analysis::Instance();

        const G4double targetMass =
            GetReactionParticleDefinition(fTargetZ, fTargetA)->GetPDGMass();
        G4ParticleDefinition *lightDefinition =
            GetReactionParticleDefinition(fLightProductZ, fLightProductA, fLightProductEx);
        G4ParticleDefinition *heavyDefinition =
            GetReactionParticleDefinition(fHeavyProductZ, fHeavyProductA, fHeavyProductEx);
        if (!lightDefinition || !heavyDefinition)
            throw std::runtime_error("Unable to create reaction product particle definitions");

        const G4double lightMass = lightDefinition->GetPDGMass();
        const G4double heavyMass = heavyDefinition->GetPDGMass();
        const G4LorentzVector projectile4 = aTrack.GetDynamicParticle()->Get4Momentum();

        TLorentzVector projectile(projectile4.px(), projectile4.py(), projectile4.pz(),
                                  projectile4.e());
        TLorentzVector target(0., 0., 0., targetMass);
        TLorentzVector cm = projectile + target;

        const G4double s = cm.M2();
        const G4double sumMass = lightMass + heavyMass;
        const G4double diffMass = lightMass - heavyMass;
        if (s < sumMass * sumMass)
            throw std::runtime_error("Configured reaction is energetically closed");

        const G4double sqrtS = std::sqrt(s);
        const G4double pCm =
            std::sqrt(std::max(0., (s - sumMass * sumMass) *
                                       (s - diffMass * diffMass))) /
            (2. * sqrtS);
        const G4double eLightCm =
            (s + lightMass * lightMass - heavyMass * heavyMass) / (2. * sqrtS);
        const G4double eHeavyCm =
            (s + heavyMass * heavyMass - lightMass * lightMass) / (2. * sqrtS);

        TVector3 lightMomentumCm(pCm, 0., 0.);
        lightMomentumCm.SetTheta(std::acos(-1. + 2. * G4UniformRand()));
        lightMomentumCm.SetPhi(G4UniformRand() * 2. * CLHEP::pi);

        TLorentzVector light4(lightMomentumCm, eLightCm);
        TLorentzVector heavy4(-lightMomentumCm, eHeavyCm);
        const TVector3 betaCm = cm.BoostVector();
        light4.Boost(betaCm);
        heavy4.Boost(betaCm);
        double ProtonPX = light4.Px();
        double ProtonPY = light4.Py();
        double ProtonPZ = light4.Pz();
        double ProtonE = (ProtonPX*ProtonPX + ProtonPY*ProtonPY + ProtonPZ*ProtonPZ)/(2*lightMass);
        double aux = atan(sqrt(ProtonPX*ProtonPX + ProtonPY*ProtonPY)/ProtonPZ) *180/CLHEP::pi;
        if(aux<0){
            double ProtonPhi = 180 + aux;
            // G4cout << ProtonPhi << G4endl;
            analysis->SetProtonPhi(ProtonPhi);
        } else{
            double ProtonPhi = aux;
            analysis->SetProtonPhi(ProtonPhi);
        }
        double ProtonTheta = acos(ProtonPX/ProtonPY) *180/CLHEP::pi;


        analysis->SetProtonPX(ProtonPX);
        analysis->SetProtonPY(ProtonPY);
        analysis->SetProtonPZ(ProtonPZ);
        analysis->SetProtonE(ProtonE);

        analysis->SetProtonTheta(ProtonTheta);


        auto *light = new G4DynamicParticle(lightDefinition,
                                            G4ThreeVector(light4.Px(), light4.Py(), light4.Pz()));
        light->SetCharge(fLightProductCharge * eplus);
        auto *heavy = new G4DynamicParticle(heavyDefinition,
                                            G4ThreeVector(heavy4.Px(), heavy4.Py(), heavy4.Pz()));
        heavy->SetCharge(fHeavyProductCharge * eplus);
        
        auto *lightTrack = new G4Track(light, aTrack.GetGlobalTime(), reactionPosition);
        auto *heavyTrack = new G4Track(heavy, aTrack.GetGlobalTime(), reactionPosition);


        aParticleChange.SetNumberOfSecondaries(3);
        aParticleChange.AddSecondary(lightTrack);
        //aParticleChange.AddSecondary(heavyTrack);
        // if(aux<0){
        ThreeBodyDecay(*heavyTrack, reactionPosition);
        // }
    }

    void ReactionProcess::ThreeBodyDecay(const G4Track &aTrack, const G4ThreeVector &reactionPosition){
        

        G4ParticleDefinition *decayDefinition =
            GetReactionParticleDefinition(fDecayHeavyZ, fDecayHeavyA, 0);
        G4ParticleDefinition *lightOneDefinition =
            GetReactionParticleDefinition(fDecayLightOneZ, fDecayLightOneA, 0);
        G4ParticleDefinition *lightTwoDefinition =
            GetReactionParticleDefinition(fDecayLightTwoZ, fDecayLightTwoZ, 0);

        const G4double decayMass = decayDefinition->GetPDGMass();
        const G4double lightOneMass = lightOneDefinition->GetPDGMass();
        const G4double lightTwoMass = lightTwoDefinition->GetPDGMass();
        const G4LorentzVector projectile4 = aTrack.GetDynamicParticle()->Get4Momentum();
        TLorentzVector projectile(projectile4.px(), projectile4.py(), projectile4.pz(),
                                  projectile4.e());
        // G4cout << (projectile4.px()*projectile4.px() + projectile4.py()*projectile4.py() + projectile4.pz()*projectile4.pz())/(2*87479.78) << G4endl;
        TLorentzVector CM = projectile;
        // This energy is correct
        // No target so I don't need to calculate the s

        const G4double s = CM.M2();
        const G4double sqrtS = std::sqrt(s);
        const G4double EkCM = fHeavyProductEx + fDecayQ;

        // This only works for two of the same decay particles ie 2n

        const G4double EMax = 4*EkCM*(lightOneMass/decayMass)/(1+4*lightOneMass/decayMass);
        const G4double ekDecayCM = EMax * G4UniformRand();
        // For single neutron emission channel VVV

        // const G4double ekDecayCM = EkCM*lightOneMass/(lightOneMass + decayMass);

        const G4double decayPCM = sqrt(2*ekDecayCM*decayMass); //Only care about tracking the heavy here
        const G4double eDecayCM = sqrt(decayPCM*decayPCM  + decayMass*decayMass);
        TVector3 decayMomentumCM(decayPCM, 0., 0.);
        decayMomentumCM.SetTheta(std::acos(-1. + 2. * G4UniformRand()));
        decayMomentumCM.SetPhi(G4UniformRand() * 2. * CLHEP::pi);
        TLorentzVector decay4(decayMomentumCM, eDecayCM);
        const TVector3 betaCm = CM.BoostVector();
        // G4cout << decay4.Pz() << G4endl;
        decay4.Boost(betaCm); //not boosting properly???
        auto *decay = new G4DynamicParticle(decayDefinition,
                                            G4ThreeVector(decay4.Px(), decay4.Py(), decay4.Pz()));
        // G4cout << decay4.Pz() << G4endl;
        decay->SetCharge(fDecayHeavyCharge * eplus);
        auto *decayTrack = new G4Track(decay, aTrack.GetGlobalTime(), reactionPosition);
        //Okay so decay is a particle actually being made, but for some reason it won't get added?
        aParticleChange.AddSecondary(decayTrack);
        // G4cout << projectile4 << G4endl;
        //Maybe??
    }
}
