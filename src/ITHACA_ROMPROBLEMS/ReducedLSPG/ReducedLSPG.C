#include "ReducedLSPG.H"

ReducedLSPG::ReducedLSPG(std::shared_ptr<ITHACAcontext> context)
: context_(std::move(context))
  // adjustTimeStep(context_->args().lookupOrDefault<bool>("adjustTimeStep", true)),
  // maxCo(context_->args().lookupOrDefault<scalar>("maxCo", 0.5)),
  // maxDeltaT(context_->args().lookupOrDefault<scalar>("maxDeltaT", 1))
{
  if (!context_)
  {
    FatalErrorInFunction
      << "Context pointer is null. Please provide a valid ITHACAcontext object."
      << exit(FatalError);
  }
}
