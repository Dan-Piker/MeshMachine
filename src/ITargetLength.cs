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
    public interface ITargetLength 
    {
        double Calculate(PlanktonMesh P, int HalfEdge);        
    }
}
