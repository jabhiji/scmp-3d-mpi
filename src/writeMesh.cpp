#include "writeMesh.h"

// this function writes mesh data to output files in XDMF format for visualization
// light data is written in XDMF format
// heavy data is written in HDF5 format

void writeMesh(const int      NX,
               const int      NY,
               const int      NZ,
               const int      time,
               const double*  rho)
{
    std::cout << "writing data to output files for t = " << time << std::endl;

    // create a 1D array of X-Y-Z coordinates for this process
    // these are "NODE CENTERED" values at the vertices of the voxels

    float *xyz = new float [NX*NY*NZ*3];

    // "natural" index of the "3D" array

    int ndx = 0; 

    // begin for loop to populate xyz

    for (int k = 0; k < NZ; k++)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int i = 0; i < NX; i++)
            {
                xyz[ndx++] = (float) i;
                xyz[ndx++] = (float) j;
                xyz[ndx++] = (float) k;
            }
        }
    }

    // CREATE FILES

    std::stringstream ss;
    ss << std::setw(6) << std::setfill('0') << time;   

    // create HDF5 files to store heavy data (mesh --> x y z coordinates of voxel vertices)
    //                                       (attribute --> fluid density, rho)

    hid_t   file_id, dataset;      // file and dataset handles
    hid_t   datatype, dataspace;   // handles
    hsize_t dimsf[1];              // dataset dimensions
    herr_t  status;

    // we will write (contiguous) 1D array data for both the mesh and for the attribute (rho)
    const int RANK = 1;

    std::string hdf5_file_with_path = "../out/hdf5/data_" + ss.str() + ".h5";
    std::string hdf5_file = "data_" + ss.str() + ".h5";
    file_id = H5Fcreate(hdf5_file_with_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // MESH DATA
    {
        // describe the size of the array and create the data space for fixed size dataset
        dimsf[0] = NX*NY*NZ*3;
        dataspace = H5Screate_simple(RANK, dimsf, NULL);

        // define datatype for the data in the file
        // we will store floating point numbers

        datatype = H5Tcopy(H5T_NATIVE_FLOAT);
        status = H5Tset_order(datatype, H5T_ORDER_LE);

        // create a new dataset within the file using defined dataspace and datatype and default dataset creation properties
    
        dataset = H5Dcreate2(file_id, "/xyz", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // write the mesh data to the dataset using default transfer properties

        status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
    }

    // NODE CENTERED DATA (rho)
    {
        // describe the size of the array and create the data space for fixed size dataset
        dimsf[0] = NX*NY*NZ;
        dataspace = H5Screate_simple(RANK, dimsf, NULL);

        // define datatype for the data in the file
        // we will store double precision numbers

        datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
        status = H5Tset_order(datatype, H5T_ORDER_LE);

        // create a new dataset within the file using defined dataspace and datatype and default dataset creation properties

        dataset = H5Dcreate2(file_id, "/rho", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // write the density data to the dataset using default transfer properties

        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho);
    }

    // release resources

    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);

    // close the HDF5 file

    H5Fclose(file_id);

    // free allocated memory

    delete [] xyz;

    // create XDMF file containing information about the mesh (light data)

    std::ofstream XDMF;
    std::string xdmf_filename = "../out/xdmf/data_" + ss.str() + ".xmf";
    XDMF.open(xdmf_filename.c_str());

    // check whether the file was opened successfully
    // if yes then continue otherwise terminate program execution

    if(XDMF.fail())
    {
        std::cout << "ERROR: could not open file for writing XDMF output!" << std::endl;
    }

    XDMF << "    <Grid Name=\"mesh " << ss.str() << "\" GridType=\"Uniform\">\n";
    XDMF << "        <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"" << NZ << " " << NY << " " << NX << "\" >\n";
    XDMF << "        </Topology>\n";
    XDMF << "        <Geometry GeometryType=\"XYZ\">\n";
    XDMF << "            <DataItem Format=\"HDF\" Precision=\" 4 \" Dimensions=\"" << NZ*NY*NX << " " << " 3 \">\n";
    XDMF << "                " << "./hdf5/" << hdf5_file << ":/xyz\n";
    XDMF << "            </DataItem>\n";
    XDMF << "        </Geometry>\n";
    XDMF << "        <Attribute Name=\"rho\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    XDMF << "            <DataItem Dimensions=\"" << NZ << " " << NY << " " << NX << "\" Precision=\" 8 \" Format=\"HDF\">\n";
    XDMF << "                " << "./hdf5/" << hdf5_file << ":/rho\n";
    XDMF << "            </DataItem>\n";
    XDMF << "        </Attribute>\n";
    XDMF << "    </Grid>\n";

    // close the file

    XDMF.close();
}
