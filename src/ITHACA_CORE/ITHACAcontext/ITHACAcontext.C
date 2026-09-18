#include "ITHACAcontext.H"

ITHACAcontext::ITHACAcontext(int argc, char *argv[])
{
    _args = autoPtr<argList>(new argList(argc, argv));

    if (!_args->checkRootCase())
    {
        Foam::FatalError.exit();
    }

    _runTime = autoPtr<Time>
    (
        new Time(Time::controlDictName, _args())

    );

    _mesh = autoPtr<fvMesh>
    (
        new fvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                _runTime->timeName(),
                _runTime(),
                IOobject::MUST_READ
            )
        )
    );
    _MRF = autoPtr<IOMRFZoneList>(new IOMRFZoneList(_mesh()));
    _pimple = autoPtr<pimpleControl>(new pimpleControl(_mesh()));
    _fvOptions = autoPtr<fv::options>(new fv::options(_mesh()));
}

ITHACAcontext::~ITHACAcontext()
{
    _fvOptions.clear();
    _pimple.clear();
    _MRF.clear();
    _mesh.clear();
    _runTime.clear();
    _args.clear();
}

argList& ITHACAcontext::args()
{
    return _args();
}

Time& ITHACAcontext::runTime()
{
    return _runTime();
}

fvMesh& ITHACAcontext::mesh()
{
    return _mesh();
}

pimpleControl& ITHACAcontext::pimple()
{
    return _pimple();
}

fv::options& ITHACAcontext::fvOptions()
{
    return _fvOptions();
}

IOMRFZoneList& ITHACAcontext::MRF()
{
    return _MRF();
}
