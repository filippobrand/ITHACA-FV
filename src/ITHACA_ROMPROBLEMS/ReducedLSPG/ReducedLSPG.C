#include "ReducedLSPG.H"

ReducedLSPG::ReducedLSPG(std::shared_ptr<ITHACAcontext> context)
: context_(std::move(context))
{
  if (!context_)
  {
    FatalErrorInFunction
      << "Context pointer is null. Please provide a valid ITHACAcontext object."
      << exit(FatalError);
  }
  ITHACAdict = new IOdictionary
  (
      IOobject
      (
          "ITHACAdict",
          runTime().system(),
          mesh(),
          IOobject::MUST_READ,
          IOobject::NO_WRITE
      )
  );

  romID_ = ITHACAdict->lookupOrDefault<word>("romID", "");
}

argList& ReducedLSPG::args() const
{
    return context_->args();
}

fvMesh& ReducedLSPG::mesh() const
{
    return context_->mesh();
}

Time& ReducedLSPG::runTime() const
{
    return context_->runTime();
}

pimpleControl& ReducedLSPG::pimple() const
{
    return context_->pimple();
}

fv::options& ReducedLSPG::fvOptions() const
{
    return context_->fvOptions();
}

IOMRFZoneList& ReducedLSPG::MRF() const
{
    return context_->MRF();
}
