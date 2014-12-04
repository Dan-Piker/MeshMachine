using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Parameters;
using Rhino.Geometry;
using Plankton;
using PlanktonGh;

namespace remesher
{
    public interface IVertexMove
    {
        List<Point3d> NewPositions(PlanktonMesh P, List<int> FeatureV, List<int> FeatureE, List<Point3d> FV, List<Curve> FC);   
    }
}
