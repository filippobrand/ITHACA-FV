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

template <typename FieldType>
PtrList<FieldType> ITHACAcontext::loadFieldsFromCaseFolder(const word& fieldName)
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    const instantList times = runTime.times();
    PtrList<FieldType> fields;

    for (const instant& time : times)
    {
        const word timeName = time.name();

        IOobject fieldIO
        (
            fieldName,
            timeName,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (!fieldIO.typeHeaderOk<FieldType>(true))
        {
            Info << "Skipping " << timeName << "/" << fieldName << nl;
            continue;
        }

        Info << "Loading field " << fieldName
             << " from time directory " << timeName << endl;

        fields.emplace_back(fieldIO, mesh);
    }

    return fields;
}

// Explicit template instantiation
template PtrList<volVectorField> ITHACAcontext::loadFieldsFromCaseFolder<volVectorField>(const word& fieldName);
template PtrList<volScalarField> ITHACAcontext::loadFieldsFromCaseFolder<volScalarField>(const word& fieldName);
template PtrList<surfaceScalarField> ITHACAcontext::loadFieldsFromCaseFolder<surfaceScalarField>(const word& fieldName);

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
