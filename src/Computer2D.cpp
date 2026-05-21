#include <memory>

#include "TGLException.h"
#include "Computer2D.h"

// inlcude header for specific computer
#include "HiCComputers.h"

// factory function
Computer2D* Computer2D::unserializeComputer2D(BufferedFile& bfile, const char *trackdb_path, int chromid1, int chromid2)
{
    Computer2D::Computer2DType type;
    bfile.read(&type, sizeof(type));

    std::unique_ptr<Computer2D> result;
    switch (type) {
    case Computer2D::CT2_AREA:
        result.reset(new AreaComputer2D(trackdb_path, chromid1, chromid2));
        break;
    case Computer2D::CT2_POTENTIAL:
        result.reset(new PotentialComputer2D(trackdb_path, chromid1, chromid2));
        break;
    case Computer2D::CT2_TECHNICAL:
        result.reset(new TechnicalComputer2D(trackdb_path, chromid1, chromid2));
        break;
    case Computer2D::CT2_TEST:
        result.reset(new TestComputer2D(trackdb_path, chromid1, chromid2));
        break;
    default:
        TGLError("Unknown computer2D type: %d", type);
    }
    // read other computer attribs; if this throws, unique_ptr releases the
    // freshly-allocated computer instead of leaking it
    result->unserialize(bfile);

    return result.release();
}

void Computer2D::serializeComputer2D(BufferedFile& bfile, Computer2D* computer)
{
    // write type
    Computer2D::Computer2DType type = computer->get_type();
	bfile.write(&type, sizeof(type));

    // write other computer attribs
    computer->serialize(bfile);
}
