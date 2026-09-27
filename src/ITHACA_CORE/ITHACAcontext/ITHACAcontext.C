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
PtrList<FieldType> ITHACAcontext::loadFieldsFromCaseFolder
(
    const word& fieldName,
    const fileName& folder
)
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    PtrList<FieldType> fields;

    const fileName searchDir = folder.empty() ? runTime.path() : runTime.path()/folder;

    // Enumerate time-like subdirectories manually, since runTime.times()
    // only knows about the case root's own time directories
    fileNameList timeDirNames = Foam::readDir(searchDir, fileName::DIRECTORY);

    DynamicList<instant> times(timeDirNames.size());
    forAll(timeDirNames, i)
    {
        if (timeDirNames[i].size() && isdigit(timeDirNames[i][0]))
        {
            times.append(instant(timeDirNames[i]));
        }
    }
    Foam::sort(times);

    for (const instant& time : times)
    {
        const word timeName = time.name();
        const fileName instance = folder.empty() ? fileName(timeName) : folder/timeName;

        IOobject fieldIO
        (
            fieldName,
            instance,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (!fieldIO.typeHeaderOk<FieldType>(true))
        {
            Info << "Skipping " << instance << "/" << fieldName << nl;
            continue;
        }

        Info << "Loading field " << fieldName
             << " from time directory " << instance << endl;
        fields.emplace_back(fieldIO, mesh);
    }

    return fields;
}

// Explicit template instantiation
template PtrList<volVectorField> ITHACAcontext::loadFieldsFromCaseFolder<volVectorField>(const word& fieldName, const fileName& folder);
template PtrList<volScalarField> ITHACAcontext::loadFieldsFromCaseFolder<volScalarField>(const word& fieldName, const fileName& folder);
template PtrList<surfaceScalarField> ITHACAcontext::loadFieldsFromCaseFolder<surfaceScalarField>(const word& fieldName, const fileName& folder);

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
