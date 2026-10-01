#ifndef NPS_PARTONS_PI0_PROJECTION_H
#define NPS_PARTONS_PI0_PROJECTION_H

#include "pi0_response_conventions.h"

// Native PARTONS v5 adapter for the diagnostic pi0 cross-section extraction.
// The including executable must be built against the installed PARTONS,
// ElementaryUtils and NumA++ libraries. No container or PDF set is selected
// here: the chosen GK19 GPD module is the one used by the installed C++ DVMP
// example, and its source does not request an LHAPDF set.

#include <ElementaryUtils/parameters/Parameter.h>
#include <ElementaryUtils/parameters/Parameters.h>
#include <partons/FundamentalPhysicalConstants.h>
#include <partons/ModuleObjectFactory.h>
#include <partons/Partons.h>
#include <partons/ServiceObjectRegistry.h>
#include <partons/beans/MesonType.h>
#include <partons/beans/PerturbativeQCDOrderType.h>
#include <partons/beans/observable/DVMP/DVMPObservableKinematic.h>
#include <partons/modules/convol_coeff_function/DVMP/DVMPCFFGK06.h>
#include <partons/modules/gpd/GPDGK19.h>
#include <partons/modules/observable/DVMP/cross_section/DVMPCrossSectionUUUMinus.h>
#include <partons/modules/process/DVMP/DVMPProcessGK06.h>
#include <partons/modules/scales/DVMP/DVMPScalesQ2Multiplier.h>
#include <partons/modules/xi_converter/DVMP/DVMPXiConverterXBToXi.h>
#include <partons/services/DVMPObservableService.h>
#include <partons/utils/type/PhysicalUnit.h>

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace nps_partons_pi0 {

struct Projection {
    bool valid = false;
    double q2 = 0.0, xb = 0.0, t = 0.0, ebeam = 0.0;
    double epsilon = 0.0;
    // FULL Hand flux, Gamma, for d4sigma/dxB dQ2 dt dphi.
    // Dividing the native nb-valued observable by this factor gives the
    // virtual-photon d2sigma/dt dphi before converting to SIMC CM units.
    double electron_flux_xbq2 = 0.0;
    // In the same raw units as the extractor's sigcm: microbarn/MeV^2.
    // The 1/(2pi) in the extractor's phi basis makes sigmaU the phi-integral
    // and sigmaLT, sigmaTT the corresponding response coefficients.
    double sigmaU = 0.0, sigmaLT = 0.0, sigmaTT = 0.0;
    // Native PARTONS electroproduction observable at three phi angles, in nb.
    std::array<double, 3> electron_nb{};
};

class Model {
public:
    Model() = default;
    Model(const Model&) = delete;
    Model& operator=(const Model&) = delete;

    ~Model() {
        // The module factory owns its objects. Release links in the same
        // order as the installed PARTONS C++ example, then close PARTONS.
        if (factory_) {
            if (gpd_) factory_->updateModulePointerReference(gpd_, 0);
            if (cff_) factory_->updateModulePointerReference(cff_, 0);
            if (xi_) factory_->updateModulePointerReference(xi_, 0);
            if (scales_) factory_->updateModulePointerReference(scales_, 0);
            if (process_) factory_->updateModulePointerReference(process_, 0);
            if (observable_) factory_->updateModulePointerReference(observable_, 0);
        }
        if (initialized_) PARTONS::Partons::getInstance()->close();
    }

    void initialize(const std::string& executable_path, int warmups, int calls) {
        if (initialized_) throw std::runtime_error("PARTONS model initialized twice.");
        if (warmups <= 0 || calls <= 0)
            throw std::runtime_error("PARTONS MC warmups and calls must be positive.");
        // PARTONS v5 locates its mandatory partons.properties beside argv[0].
        // The pipeline creates that small configuration beside the xsec binary.
        std::vector<char> argv0(executable_path.begin(), executable_path.end());
        argv0.push_back('\0');
        char* partons_argv[] = {argv0.data()};
        PARTONS::Partons* app = PARTONS::Partons::getInstance();
        app->init(1, partons_argv);
        initialized_ = true;

        factory_ = app->getModuleObjectFactory();
        service_ = app->getServiceObjectRegistry()->getDVMPObservableService();
        gpd_ = factory_->newGPDModule(PARTONS::GPDGK19::classId);
        cff_ = factory_->newDVMPConvolCoeffFunctionModule(PARTONS::DVMPCFFGK06::classId);

        ElemUtils::Parameters cff_parameters;
        cff_parameters.add(ElemUtils::Parameter(
            PARTONS::PerturbativeQCDOrderType::PARAMETER_NAME_PERTURBATIVE_QCD_ORDER_TYPE,
            PARTONS::PerturbativeQCDOrderType::LO));
        cff_parameters.add(ElemUtils::Parameter(
            PARTONS::DVMPCFFGK06::PARAMETER_NAME_DVMPCFFGK06_MC_NWARMUP, warmups));
        cff_parameters.add(ElemUtils::Parameter(
            PARTONS::DVMPCFFGK06::PARAMETER_NAME_DVMPCFFGK06_MC_NCALLS, calls));
        cff_parameters.add(ElemUtils::Parameter(
            PARTONS::DVMPCFFGK06::PARAMETER_NAME_DVMPCFFGK06_MC_CHI2LIMIT, 0.8));
        cff_->configure(cff_parameters);

        xi_ = factory_->newDVMPXiConverterModule(PARTONS::DVMPXiConverterXBToXi::classId);
        scales_ = factory_->newDVMPScalesModule(PARTONS::DVMPScalesQ2Multiplier::classId);
        scales_->configure(ElemUtils::Parameter(
            PARTONS::DVMPScalesQ2Multiplier::PARAMETER_NAME_LAMBDA, 1.0));
        process_ = factory_->newDVMPProcessModule(PARTONS::DVMPProcessGK06::classId);
        observable_ = factory_->newDVMPObservable(PARTONS::DVMPCrossSectionUUUMinus::classId);
        observable_->setProcessModule(process_);
        process_->setScaleModule(scales_);
        process_->setXiConverterModule(xi_);
        process_->setConvolCoeffFunctionModule(cff_);
        cff_->setGPDModule(gpd_);
    }

    Projection predict(double q2, double xb, double t, double ebeam) const {
        if (!initialized_ || !service_ || !observable_)
            throw std::runtime_error("PARTONS model is not initialized.");
        Projection out;
        out.q2 = q2; out.xb = xb; out.t = t; out.ebeam = ebeam;
        if (!(std::isfinite(q2) && q2 > 0.0 && std::isfinite(xb) && xb > 0.0 && xb < 1.0 &&
              std::isfinite(t) && t < 0.0 && std::isfinite(ebeam) && ebeam > 0.0))
            return out;

        const double mp = PARTONS::Constant::PROTON_MASS;
        const double pi = PARTONS::Constant::PI;
        const double y = q2 / (2.0 * mp * xb * ebeam);
        const double gamma = 2.0 * xb * mp / std::sqrt(q2);
        const double target_mass_term = std::pow(y * gamma / 2.0, 2.0);
        out.epsilon = (1.0 - y - target_mass_term) /
                      (1.0 - y + y * y / 2.0 + target_mass_term);
        if (!(std::isfinite(out.epsilon) && out.epsilon > 0.0 && out.epsilon < 1.0))
            return out;

        // Installed DVMPProcessGK06::CrossSection starts with
        // alpha*(W2-M2)/(16*pi2*E2*M2*Q2*(1-eps)), then applies Q2/xB2.
        // This is Gamma/(2*pi), multiplying the phi-integrated response
        // bracket, NOT the full Hand flux. Its final target-azimuth 1/(2*pi)
        // is canceled inside DVMPCrossSectionUUUMinus::computeObservable.
        // Use Gamma here because the projection below explicitly restores
        // the hadron-azimuth 2*pi when extracting U, LT, and TT.
        out.electron_flux_xbq2 = nps_pi0_conventions::hand_flux_xbq2(
            q2, xb, ebeam, out.epsilon, mp, PARTONS::Constant::FINE_STRUCTURE_CONSTANT);
        if (!(std::isfinite(out.electron_flux_xbq2) && out.electron_flux_xbq2 > 0.0))
            return out;

        // The UUMinus observable averages beam helicities. Its phi shape has
        // only 1, cos(phi), and cos(2phi) for an unpolarized target. Three
        // angles therefore determine all coefficients. PARTONS caches the
        // CFF for identical xB, Q2, t and scales, so changing only phi does
        // not rerun the MC convolution. Numerical MC uncertainty is not
        // exposed by this observable and must not be represented as zero.
        const std::array<double, 3> phis{0.0, pi / 2.0, pi};
        for (size_t i = 0; i < phis.size(); ++i) {
            PARTONS::DVMPObservableKinematic kin(
                xb, t, q2, ebeam, phis[i], PARTONS::MesonType::PI0);
            const auto result = service_->computeSingleKinematic(kin, observable_);
            out.electron_nb[i] = result.getValue().makeSameUnitAs(PARTONS::PhysicalUnit::NB).getValue();
            if (!(std::isfinite(out.electron_nb[i]) && out.electron_nb[i] >= 0.0))
                return out;
        }

        // Keep the flux, phi convention, and unit conversion independently
        // testable without invoking the stochastic GK convolution.
        const auto responses = nps_pi0_conventions::project_electron_observables(
            out.electron_nb, out.electron_flux_xbq2, out.epsilon);
        out.sigmaU = responses.U;
        out.sigmaLT = responses.LT;
        out.sigmaTT = responses.TT;
        // DVMPProcessGK06::CrossSection returns an exact zero when a point
        // lies between its GK approximation to t_min and the exact forward
        // limit ("|t| smaller than that used by GK"). That is a model-domain
        // rejection, not a prediction of vanishing pi0 production. Require
        // a positive unpolarized response so the extractor marks it absent.
        out.valid = responses.valid;
        return out;
    }

private:
    bool initialized_ = false;
    PARTONS::ModuleObjectFactory* factory_ = nullptr;
    PARTONS::DVMPObservableService* service_ = nullptr;
    PARTONS::GPDModule* gpd_ = nullptr;
    PARTONS::DVMPConvolCoeffFunctionModule* cff_ = nullptr;
    PARTONS::DVMPXiConverterModule* xi_ = nullptr;
    PARTONS::DVMPScalesModule* scales_ = nullptr;
    PARTONS::DVMPProcessModule* process_ = nullptr;
    PARTONS::DVMPObservable* observable_ = nullptr;
};

} // namespace nps_partons_pi0

#endif
