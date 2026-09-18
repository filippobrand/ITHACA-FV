#include "ReducedLSPG.H"

ReducedLSPG::ReducedLSPG(int argc, char* argv[])
{
  _args = autoPtr<argList>(
                new argList(argc, argv));

  if (!_args->checkRootCase())
  {
      Foam::FatalError.exit();
  }
  argList& args = _args();
  #include "createTime.H"
  #include "createMesh.H"
  _MRF = autoPtr<IOMRFZoneList>(
           new IOMRFZoneList(mesh));
  _pimple = autoPtr<pimpleControl>(
                new pimpleControl(
                    mesh));
  #include "createFvOptions.H"
}
