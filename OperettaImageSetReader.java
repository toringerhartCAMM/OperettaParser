/*
 * #%L
 * OME Bio-Formats package for reading and converting biological file formats.
 * %%
 * Copyright (C) 2005 - 2016 Open Microscopy Environment:
 *   - Board of Regents of the University of Wisconsin-Madison
 *   - Glencoe Software, Inc.
 *   - University of Dundee
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 2 of the 
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public 
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 * #L%
 */

package loci.formats.in;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.HashMap;

import loci.common.DataTools;
import loci.common.Location;
import loci.common.RandomAccessInputStream;
import loci.common.xml.BaseHandler;
import loci.common.xml.XMLTools;

import loci.formats.CoreMetadata;
import loci.formats.FormatException;
import loci.formats.FormatReader;
import loci.formats.FormatTools;
import loci.formats.MetadataTools;
import loci.formats.meta.MetadataStore;
import loci.formats.tiff.IFD;
import loci.formats.tiff.TiffParser;
import loci.formats.tiff.TiffIFDEntry;

import ome.units.UNITS;
import ome.units.quantity.Length;
import ome.units.quantity.Time;
import ome.units.quantity.Temperature;

import ome.xml.model.primitives.PercentFraction;
import ome.xml.model.primitives.NonNegativeInteger;
import ome.xml.model.primitives.PositiveFloat;
import ome.xml.model.primitives.PositiveInteger;

import org.xml.sax.Attributes;


public class OperettaImageSetReader extends FormatReader {

    // -- Constants --
    private static final int XML_TAG  = 65500;
    private static final int USER_TAG = 315;

    // -- Fields --
    private MinimalTiffReader reader;

    /** List of all files to open */
    private List<String> allFiles;
    private String fileLocation;
    private String[][] orderedFiles;
    private ImageInfo[][] allImages;

    // -- Constructor --
    /** Constructs a new Operetta reader. */
    public OperettaImageSetReader() {
        super("PerkinElmer Operetta ImageSet", new String[] { "tif", "tiff" });
        domains = new String[] { FormatTools.HCS_DOMAIN };
        suffixNecessary = true;
        suffixSufficient = false;
        hasCompanionFiles = true;
        datasetDescription = "A directory containing relocated .tif/.tiff files from the Operetta.";
    }


    // -- IFormatReader API methods --
    /* @see loci.formats.IFormatReader#getRequiredDirectories(String[]) */
    @Override
    public int getRequiredDirectories(String[] files) throws FormatException, IOException {
        return 1;
    }


    /* @see loci.formats.IFormatReader#isSingleFile(String) */
    @Override
    public boolean isSingleFile(String id) throws FormatException, IOException {
        return false;
    }


    /* @see loci.formats.IFormatReader#fileGroupOption(String) */
    @Override
    public int fileGroupOption(String id) throws FormatException, IOException {
        return FormatTools.MUST_GROUP;
    }


    /* @see loci.formats.IFormatReader#isThisType(String, boolean) */
    @Override
    public boolean isThisType(String name, boolean open) {
        String localName = new Location(name).getName();
        if (localName.endsWith(".tif") || localName.endsWith(".tiff")) {
            return true;
        }
        return false;
    }


    /* @see loci.formats.IFormatReader#isThisType(RandomAccessInputStream) */
    @Override
    public boolean isThisType(RandomAccessInputStream stream) throws IOException {
        TiffParser p = new TiffParser(stream);
        IFD ifd = p.getFirstIFD();
        if (ifd == null) return false;

        Object s = ifd.getIFDValue(XML_TAG);
        if (s == null) return false;
        
        String xml = s instanceof String[] ? ((String[]) s)[0] : s.toString();
        return xml.indexOf("Operetta") < 1024;
    }


    /* @see loci.formats.IFormatReader#close(boolean) */
    @Override
    public void close(boolean fileOnly) throws IOException {
        super.close(fileOnly);
        if (!fileOnly) {
            if (reader != null) {
                reader.close();
            }
            reader = null;
            allFiles = null;
            orderedFiles = null;
            allImages = null;
            fileLocation = null;
        }
    }


    // @see loci.formats.IFormatReader#openBytes(int, byte[], int, int, int, int)
    @Override
    public byte[] openBytes(int no, byte[] buf, int x, int y, int w, int h)
        throws FormatException, IOException
    {
        FormatTools.checkPlaneParameters(this, no, buf.length, x, y, w, h);
        
        if (getSeries() < orderedFiles.length && no < orderedFiles[getSeries()].length) {
            String filename = orderedFiles[getSeries()][no];
            
            if (new Location(filename).exists()) {
                if (reader == null) {
                    reader = new MinimalTiffReader();
                }

                reader.setId(filename);
                reader.openBytes(0, buf, x, y, w, h);
                reader.close();
            }
        }
        
        return buf;
    }


    /* @see loci.formats.IFormatReader#getSeriesUsedFiles(boolean) */
    @Override
    public String[] getSeriesUsedFiles(boolean noPixels) {
        FormatTools.assertId(currentId, true, 1);
        
        ArrayList<String> files = new ArrayList<String>();
        files.add(currentId);
        for (String filename : orderedFiles[getSeries()]) {
            files.add(filename);
        }
        
        return files.toArray(new String[files.size()]);
    }


    public OperettaXMLHandler parseXmlFromTiffHeader (String file)
        throws IOException
    {
        Location f;
        if (fileLocation == null) {
            f = new Location(file);
        } else {
            f = new Location(fileLocation, file);
        }
        
        RandomAccessInputStream s = new RandomAccessInputStream(f.getPath(), 32);
        TiffParser parser = new TiffParser(s);
    
        parser.setDoCaching(false);
        
        IFD firstIFD = parser.getFirstIFD();
        if (firstIFD == null) {
            throw new IOException("Could not read IFD in " + file + "\n");
        }

        // Read the Perkin-Elmer/Operetta tag
        Object ifd_entry = firstIFD.getIFDValue(XML_TAG);
        if (ifd_entry == null) {
            throw new IOException("Could not find Perkin-Elmer Operetta tag 65500 in " + file + "\n");
        }

        // Read the IDF entry
        TiffIFDEntry entry = TiffIFDEntry.class.cast(ifd_entry);
        int size = entry.getValueCount();
        byte[] byte_xml = new byte[size];

        // Seek to the given offset and read the data
        s.seek(entry.getValueOffset());
        s.read(byte_xml);

        // Parse the XML data
        OperettaXMLHandler tiffhandler = new OperettaXMLHandler();
        XMLTools.parseXML(byte_xml, tiffhandler);
        s.close();
        return tiffhandler; 
    }
    
    
    /* @see loci.formats.FormatReader#initFile(String) */
    @Override
    protected void initFile (String id) 
        throws FormatException, IOException 
    {
        if (!checkSuffix(id, "tif") && !checkSuffix(id, "tiff")) {
            throw new FormatException("Import started on a non-tif/tiff file.");
        }
        
        super.initFile(id);
        
        allFiles = new ArrayList<String>();
        OperettaXMLHandler tiffhandler = parseXmlFromTiffHeader(id);
        ExposureInfo exposureInfo = tiffhandler.getExposureInfo();
        
        ArrayList<StackPlaneEntry> planes = tiffhandler.getPlanes();
        ArrayList<WellEntry> wells = tiffhandler.getWells();
        ArrayList<FieldEntry> fields = tiffhandler.getFieldLocations();
        ArrayList<ExposureRecordEntry> exposures = tiffhandler.getExposures();
        ArrayList<Repetition> repetitions = tiffhandler.getKineticSequences();
        ArrayList<ObjectiveInfo> objectives = tiffhandler.getObjectives();
        ArrayList<LightSourceInfo> lightSources = tiffhandler.getLightSources();
        ArrayList<FilterChangerInfo> filterChangers = tiffhandler.getFilterChangers();
        String username = null;
        
        // Check for null/empty
        if ((planes.isEmpty()) || (planes == null)) {
            throw new IOException("StackPlaneEntry: planes is empty or null");
        }
        if ((wells.isEmpty()) || (wells == null)) {
            throw new IOException("WellEntry: wells is empty or null");
        }
        if ((fields.isEmpty()) || (fields == null)) {
            throw new IOException("FieldEntry: fields is empty or null");
        }
        if ((exposures.isEmpty()) || (exposures == null)) {
            throw new IOException("ExposureEntry: exposures is empty or null");
        }
        if ((repetitions.isEmpty()) || (repetitions == null)) {
            throw new IOException("Repetition: repetitions is empty or null");
        }
        
        // Similar to OperettaReader.java
        ArrayList<Integer> uniquerows = new ArrayList<Integer>();
        ArrayList<Integer> uniquecolumns = new ArrayList<Integer>();
        ArrayList<Integer> uniquefields = new ArrayList<Integer>();
        ArrayList<Integer> uniquezs = new ArrayList<Integer>();
        ArrayList<Integer> uniquets = new ArrayList<Integer>();
        ArrayList<Integer> uniquecs = new ArrayList<Integer>();
        
        // Now populate...
        for (WellEntry well: wells) {
            if (!uniquerows.contains(well.row)) {
                uniquerows.add(well.row);
            }
            if (!uniquecolumns.contains(well.column)) {
                uniquecolumns.add(well.column);
            }
        }
        for (Integer i = 1; i <= fields.size(); i++) {
            if (!uniquefields.contains(i)) {
                uniquefields.add(i);
            }
        }
        for (Integer i = 1; i <= exposures.size(); i++) {
            if (!uniquecs.contains(i)) {
                uniquecs.add(i);
            }
        }
        for (Integer i = 1; i <= planes.size(); i++) {
            if (!uniquezs.contains(i)) {
                uniquezs.add(i);
            }
        }
        for (Integer i = 1; i <= repetitions.size(); i++) {
            if (!uniquets.contains(i)) {
                uniquets.add(i);
            }
        }
      
        // Info arrays, then sort
        Integer[] rows = uniquerows.toArray(new Integer[uniquerows.size()]);
        Integer[] columns = uniquecolumns.toArray(new Integer[uniquecolumns.size()]);
        Integer[] allfields = uniquefields.toArray(new Integer[uniquefields.size()]);
        Integer[] zs = uniquezs.toArray(new Integer[uniquezs.size()]);
        Integer[] ts = uniquets.toArray(new Integer[uniquets.size()]);
        Integer[] cs = uniquecs.toArray(new Integer[uniquecs.size()]);
        
        Arrays.sort(rows);
        Arrays.sort(columns);
        Arrays.sort(allfields);
        Arrays.sort(zs);
        Arrays.sort(ts);
        Arrays.sort(cs);
        
        // seriesCount - total number of fields (in the entire plate)
        // planesCount - total number of planes per field
        // total number of images = seriesCount * planesCount
        int seriesCount = rows.length * columns.length * allfields.length;
        int planesCount = zs.length * cs.length * ts.length;
        
        // Get a list of all the files
        if (allFiles.isEmpty()) {
            Location parent = new Location(id).getAbsoluteFile().getParentFile();
            fileLocation = parent.toString();
            String[] fileList = parent.list(false);
            
            for (String f: fileList) {
                if ((checkSuffix(f, "tif") || checkSuffix(f, "tiff")) && !allFiles.contains(f)) {
                    allFiles.add(f);
                }
            }
        }

        // Init some variables before we read xml from every file
        Temperature[][] temperature = new Temperature[seriesCount][planesCount];
        PercentFraction[][] co2concentration = new PercentFraction[seriesCount][planesCount];
        
        /* 
         * imginfo - key:value pairs correspoding to index:imageInfo
         * 
         * where, index     - index of the image in the file list allFiles
         *        imageInfo - the metadata that describes row, column, field, 
         *                    z, t, c, of the image.
         * 
         * The metadata from imageInfo is used to determine the image's unique 
         * "imageIndex" in the dataset. 
         * 
         */
        Map<Integer,ImageInfo> imginfo = new HashMap<Integer,ImageInfo>();
        OperettaXMLHandler tiffxml;
        for (int i = 0; i < allFiles.size(); i++) {
            tiffxml = parseXmlFromTiffHeader(allFiles.get(i));
            imginfo.put(i, tiffxml.getImageInfo());
        }
        
        // Init variables
        allImages = new ImageInfo[seriesCount][planesCount];
        orderedFiles = new String[seriesCount][planesCount];
        
        // Put everything in the expected order
        int nextSeries = 0;
        for (int row = 0; row < rows.length; row++) {
            for (int col = 0; col < columns.length; col++) {
                for (int field = 0; field < allfields.length; field++) {
                    
                    int nextPlane = 0;
                    
                    for (int t = 0; t < ts.length; t++) {
                        for (int z = 0; z < zs.length; z++) {
                            for (int c = 0; c < cs.length; c++) {
                                
                                Boolean found = false;
                                
                                for (Integer index: imginfo.keySet()) {
                                    ImageInfo info = imginfo.get(index);
                                    
                                    if (info.row.equals(rows[row]) &&
                                        info.column.equals(columns[col]) &&
                                        info.field.equals(allfields[field]) &&
                                        info.plane.equals(zs[z]) &&
                                        info.record.equals(cs[c]) &&
                                        info.kinetic.equals(ts[t])) 
                                        {
                                            allImages[nextSeries][nextPlane] = info;
                                            orderedFiles[nextSeries][nextPlane] = fileLocation + "/" + allFiles.get(index);
                                            found = true;
                                            break;
                                        }
                                }
                                nextPlane++;
                            }
                        }
                    }
                    nextSeries++;
                }
            }
        }     

        Integer objectiveMagnification = null;
        for (ObjectiveInfo objective: objectives) {
            if (objective.ID.equals(exposureInfo.objectiveName)) {
                objectiveMagnification = objective.magnification;
                break;
            }
        }
        if (objectiveMagnification == null) {
            throw new IOException("No objective magnification found.");
        }
        
        core.clear();
        reader = new MinimalTiffReader();
        
        // Populate coremetadata 
        for (int i = 0; i < seriesCount; i++)
        {
            CoreMetadata ms = new CoreMetadata();
            core.add(ms);
            
            ms.sizeX = tiffhandler.getCameraXPixels();
            ms.sizeY = tiffhandler.getCameraYPixels();
            ms.sizeZ = uniquezs.size();
            ms.sizeC = uniquecs.size();
            ms.sizeT = uniquets.size();
            ms.dimensionOrder = "XYCZT";
            ms.rgb = false;
            ms.imageCount = getSizeZ() * getSizeC() * getSizeT();
            
            for (int j = 0; j < planesCount; j++)
            {
                RandomAccessInputStream s;
                Location f;
                
                f = new Location(orderedFiles[i][j]);

                try {
                    s = new RandomAccessInputStream(f.getPath(), 32);
                } catch (FileNotFoundException e) {
                    throw new FileNotFoundException("File=" + orderedFiles[i][j] + ", Bad index: i=" + i + ", j=" + j + "\n");
                }
                
                TiffParser parser = new TiffParser(s);
                parser.setDoCaching(false);

                IFD firstIFD = parser.getFirstIFD();
                if (firstIFD == null) continue;
                
                Object ifd_entry = firstIFD.getIFDValue(XML_TAG);
                if (ifd_entry == null) continue;
                
                TiffIFDEntry entry = TiffIFDEntry.class.cast(ifd_entry);
                int size = entry.getValueCount();
                byte[] byte_xml = new byte[size];
                
                s.seek(entry.getValueOffset());
                s.read(byte_xml);

                OperettaXMLHandler imgxml = new OperettaXMLHandler();
                XMLTools.parseXML(byte_xml, imgxml);
                
                temperature[i][j] = imgxml.getCurrentTemperature();
                co2concentration[i][j] = imgxml.getCurrentCO2();
                
                if (j == 0) {
                    ms.littleEndian = firstIFD.isLittleEndian();
                    ms.pixelType = firstIFD.getPixelType();
                }
                
                // User name
                Object ifd_user_entry = firstIFD.getIFDValue(USER_TAG);
                if (ifd_entry == null) continue;
                
                TiffIFDEntry userentry = TiffIFDEntry.class.cast(ifd_user_entry);
                size = userentry.getValueCount();
                byte[] username_bytes = new byte[size];
                
                s.seek(userentry.getValueOffset());
                s.read(username_bytes);
                username = new String(username_bytes);
                
                s.close();
                
            }
        }
        
        // populate the MetadataStore
        MetadataStore store = makeFilterMetadata();
        MetadataTools.populatePixels(store, this, true);
        
        // Init IDs
        String instrumentID   = MetadataTools.createLSID("Instrument", 0);
        String plateAcquID    = MetadataTools.createLSID("PlateAcquisition", 0, 0);
        String plateID        = MetadataTools.createLSID("Plate", 0);
        String experimenterID = MetadataTools.createLSID("Experimenter", 0);
        
        // Populate metadatastore
        store.setMicroscopeManufacturer("Perkin-Elmer", 0);
        store.setMicroscopeModel("Operetta", 0);
        store.setMicroscopeSerialNumber(tiffhandler.getSerial(), 0);
        store.setInstrumentID(instrumentID, 0);
        store.setPlateID(plateID, 0);
        store.setPlateAcquisitionID(plateAcquID, 0, 0);
        store.setExperimenterID(experimenterID, 0);
        store.setPlateRows(new PositiveInteger(tiffhandler.getPlateEndRow()), 0);
        store.setPlateColumns(new PositiveInteger(tiffhandler.getPlateEndColumn()), 0);
        store.setPlateExternalIdentifier(tiffhandler.getPlateID(), 0);
        store.setPlateAcquisitionName(tiffhandler.getMeasurementID(), 0, 0);
        store.setPlateDescription(tiffhandler.getPlateDescription(), 0);
        
        PositiveInteger fieldCount = FormatTools.getMaxFieldCount(allfields.length);
        if (fieldCount != null) {
            store.setPlateAcquisitionMaximumFieldCount(fieldCount, 0, 0);
        }
        
        for (int row = 0; row < rows.length; row++) {
            for (int col = 0; col < columns.length; col++) {
                
                int well = row * columns.length + col;
                
                store.setWellID(MetadataTools.createLSID("Well", 0, well), 0, well);
                store.setWellRow(new NonNegativeInteger(rows[row]), 0, well);
                store.setWellColumn(new NonNegativeInteger(columns[col]), 0, well);
                
                for (int field = 0; field < allfields.length; field++) {
                    
                    int imageIndex = well * allfields.length + field;
                    
                    String pixelsID = MetadataTools.createLSID("Pixels", imageIndex);
                    String imageID  = MetadataTools.createLSID("Image", imageIndex);
                    String wellSampleID = MetadataTools.createLSID("WellSample", 0, well, field);
                    String name = "Well " + (well + 1) + ", Field " + (field + 1);
                                        
                    store.setWellSampleID(wellSampleID, 0, well, field);
                    store.setWellSampleIndex(new NonNegativeInteger(imageIndex), 0, well, field);
                    store.setWellSampleImageRef(imageID, 0, well, field);
                    store.setPixelsID(pixelsID, imageIndex);
                    store.setImageID(imageID, imageIndex); 
                    store.setImageInstrumentRef(instrumentID, imageIndex);
                    store.setImageName(name, imageIndex); 
                    store.setPlateAcquisitionWellSampleRef(wellSampleID, 0, 0, imageIndex); 
                    
                    for (int c = 0; c < getSizeC(); c++) {
                        ExposureRecordEntry exposure = exposures.get(allImages[imageIndex][c].record - 1);
                        
                        store.setChannelName(exposure.channelName, imageIndex, c);
                        
                        Integer excitationWavelength = null;
                        for (LightSourceInfo lightSource: lightSources) {
                            if (lightSource.ID.equals(exposure.lightSourceID)) {
                                for (LightSourceOutputInfo lightSourceOutput: lightSource.outputs) {
                                    if (lightSourceOutput.ID.equals(exposure.lightSourceOutputName)) {
                                        if (lightSourceOutput.wavelength > 0) {
                                            store.setChannelExcitationWavelength(new Length(lightSourceOutput.wavelength, UNITS.NM), imageIndex, c);
                                        }
                                    }
                                }
                            }
                        }
                        
                        Integer emissionWavelength = null;
                        for (FilterChangerInfo filterChanger: filterChangers) {
                            if (filterChanger.ID.equals(exposure.emissionFilterID)) {
                                for (FilterInfo filter: filterChanger.filters) {
                                    if (filter.ID.equals(exposure.emissionFilterName)) {
                                        if (filter.wavelength > 0) {
                                            store.setChannelEmissionWavelength(new Length(filter.wavelength, UNITS.NM), imageIndex, c);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    Double avgTemp = 0.0;
                    Double avgCO2  = 0.0;
                    for (int p = 0; p < getImageCount(); p++) {
                        FieldEntry thisfield = fields.get(allImages[imageIndex][p].field - 1);
                        StackPlaneEntry plane = planes.get(allImages[imageIndex][p].plane - 1);
                        ExposureRecordEntry exposure = exposures.get(allImages[imageIndex][p].record - 1);
                        
                        store.setPlanePositionX(new Length(thisfield.x * 1000000, UNITS.REFERENCEFRAME), imageIndex, p);
                        store.setPlanePositionY(new Length(thisfield.y * 1000000, UNITS.REFERENCEFRAME), imageIndex, p);
                        store.setPlanePositionZ(new Length(plane.z * 1000000, UNITS.REFERENCEFRAME), imageIndex, p);
                        store.setPlaneExposureTime(new Time(exposure.exposureTime, UNITS.SECOND), imageIndex, p); 
                        
                        avgTemp += temperature[imageIndex][p].value().doubleValue();
                        avgCO2 += co2concentration[imageIndex][p].getValue();
                    }
                    // Set average temp and co2
                    avgTemp /= planesCount;
                    avgCO2 /= planesCount;
                    store.setImagingEnvironmentTemperature(new Temperature(avgTemp, UNITS.CELSIUS), imageIndex);
                    store.setImagingEnvironmentCO2Percent(new PercentFraction(avgCO2.floatValue()), imageIndex);
                    
                    ExposureRecordEntry exposure = exposures.get(allImages[imageIndex][0].record - 1);
                    double physicalsizex = (tiffhandler.getCameraXPixelSize() * exposure.binningx) / (objectiveMagnification * tiffhandler.getMagnification());
                    double physicalsizey = (tiffhandler.getCameraYPixelSize() * exposure.binningy) / (objectiveMagnification * tiffhandler.getMagnification());
                    store.setPixelsPhysicalSizeX(FormatTools.getPhysicalSizeX(physicalsizex * 1000000), imageIndex);
                    store.setPixelsPhysicalSizeY(FormatTools.getPhysicalSizeY(physicalsizey * 1000000), imageIndex);
                }
            }
        }
    }


    // -- Helper classes --
    class OperettaXMLHandler extends BaseHandler {
        
        private PlateInfo plateInfo;
        private ImageInfo imageInfo;
        private SublayoutInfo sublayoutInfo;
        private FieldEntry fieldEntry;
        private WellEntry wellEntry;
        private StackPlaneEntry planeEntry;
        private ExposureInfo exposureInfo;
        private CameraInfo cameraInfo;
        private ExposureRecordEntry exposureRecordEntry;
        private SlowKineticInfo slowKineticInfo;
        private ObjectiveInfo objectiveInfo;
        private LightSourceInfo lightSourceInfo;
        private LightSourceOutputInfo lightSourceOutputInfo;
        private FilterChangerInfo filterChangerInfo;
        private FilterInfo filterInfo;
        private Repetition repetition;
        
        private String currentName;
        private String displayName;
        private String plateID;
        private String measurementTime;
        private String plateName;
        private String plateDescription;
        private String measurementID;
        private String serial;
        private Double magnification;

        private Temperature currentTemperature;
        private Temperature targetTemperature;
        private PercentFraction currentCO2;
        private PercentFraction targetCO2;

        private ArrayList<StackPlaneEntry> planes = new ArrayList<StackPlaneEntry>();
        private ArrayList<WellEntry> wells = new ArrayList<WellEntry>();
        private ArrayList<ObjectiveInfo> objectives = new ArrayList<ObjectiveInfo>();
        private ArrayList<FilterChangerInfo> filterChangers = new ArrayList<FilterChangerInfo>();
        private ArrayList<LightSourceInfo> lightSources = new ArrayList<LightSourceInfo>();
        
        private StringBuffer currentValue = new StringBuffer();
        
        private Boolean plateTypeSection = false;
        private Boolean sublayoutSection = false;
        private Boolean imageInfoSection = false;
        private Boolean stackPlanesSection = false;
        private Boolean exposureSection = false;
        private Boolean cameraSection = false;
        private Boolean slowKineticSection = false;
        private Boolean objectiveSection = false;
        private Boolean lightSourceSection = false;
        private Boolean lightSourceOutputSection = false;
        private Boolean filterChangerSection = false;
        private Boolean filterSection = false;
        
        public int getImageWellRow () { 
            return imageInfo.row;
        }
        public int getImageWellColumn () { 
            return imageInfo.column;
        }
        public int getImageFieldID () { 
            return imageInfo.field;
        }
        public int getImagePlaneID () { 
            return imageInfo.plane;
        }
        public int getImageSlowKineticID () { 
            return imageInfo.kinetic;
        }
        public int getImageExposureID () { 
            return imageInfo.exposure;
        }
        public int getImageRecordID () { 
            return imageInfo.record;
        }
        public int getImageChannelID () { 
            return imageInfo.channel;
        }
        public int getImageFlimID () { 
            return imageInfo.flim;
        }
        public String getImageAbsTime () { 
            return imageInfo.absTime;
        }
        public Temperature getCurrentTemperature () {
            return currentTemperature;
        }
        public PercentFraction getCurrentCO2 () {
            return currentCO2;
        }
        public String getSerial () {
            return serial;
        }
        public Double getMagnification () {
            return magnification;
        }
        public String getMeasurementID () {
            return measurementID;
        }
        public String getLens () {
            return sublayoutInfo.lensName;
        }
        public double getOverlapX () {
            return sublayoutInfo.overlapX;
        }
        public double getOverlapY () {
            return sublayoutInfo.overlapY;
        }
        public ArrayList<FieldEntry> getFieldLocations () {
            return sublayoutInfo.fieldEntries;
        }
        public int getFieldCount () {
            return sublayoutInfo.fieldEntries.size();
        }
        public String getPlateID () {
            return plateInfo.plateId;
        }
        public String getPlateDescription () {
            return plateInfo.plateDescription;
        }
        public String getPlateWellFormTop () {
            return plateInfo.wellFormTop;
        }
        public String getPlateWellFormBottom () {
            return plateInfo.wellFormBottom;
        }
        public int getPlateStartColumn () {
            return plateInfo.start_column;
        }
        public int getPlateEndColumn () {
            return plateInfo.end_column;
        }
        public int getPlateStartRow () {
            return plateInfo.start_row;
        }
        public int getPlateEndRow () {
            return plateInfo.end_row;
        }
        public Double getPlateStartX () {
            return plateInfo.start_x;
        }
        public Double getPlateStartY () {
            return plateInfo.start_y;
        }
        public Double getPlateEndX () {
            return plateInfo.end_x;
        }
        public Double getPlateEndY () {
            return plateInfo.end_y;
        }
        public Double getFootPrintWidthX () {
            return plateInfo.footprint_width_x;
        }
        public Double getFootPrintWidthy () {
            return plateInfo.footprint_width_y;
        }
        public ArrayList<WellEntry> getWells () {
            return wells;
        }
        public ArrayList<StackPlaneEntry> getPlanes () {
            return planes;
        }
        public ArrayList<ExposureRecordEntry> getExposures () {
            return exposureInfo.exposureRecordEntries;
        } 
        public String getObjectiveName () {
            return exposureInfo.objectiveName;
        }
        public Integer getCameraXPixels () {
            return cameraInfo.xpixels;
        }
        public Integer getCameraYPixels () {
            return cameraInfo.ypixels;
        }
        public Double getCameraXPixelSize () {
            return cameraInfo.xpixelsize;
        }
        public Double getCameraYPixelSize () {
            return cameraInfo.ypixelsize;
        }
        public ExposureInfo getExposureInfo () {
            return exposureInfo;
        }
        public ArrayList<Repetition> getKineticSequences () {
            return slowKineticInfo.timepoints;
        }
        public ImageInfo getImageInfo () {
            return imageInfo;
        }
        public ArrayList<ObjectiveInfo> getObjectives () {
            return objectives;
        }
        public FilterChangerInfo getFilterChangerInfo () {
            return filterChangerInfo;
        }
        public ArrayList<LightSourceOutputInfo> getLightSourceOutputInfo () {
            return lightSourceInfo.outputs;
        }
        public LightSourceInfo getLightSourceInfo () {
            return lightSourceInfo;
        }
        public ArrayList<LightSourceInfo> getLightSources () {
            return lightSources;
        }
        public ArrayList<FilterChangerInfo> getFilterChangers () {
            return filterChangers;
        }
        
        @Override
        public void characters(char[] ch, int start, int length) {
            String value = new String(ch, start, length);
            currentValue.append(value);
        }
        
        @Override
        public void startElement(String uri, String localName, String qName, Attributes attributes)
        {
            currentValue.setLength(0);
            currentName = qName;

            if (qName.equals("Image") && attributes.getValue("id") == null && imageInfoSection == false) 
            {
                imageInfo = new ImageInfo();
                imageInfoSection = true;
            } 
            else if (qName.equals("PlateType") && attributes.getValue("id") == null && plateTypeSection == false) 
            {
                plateInfo = new PlateInfo();
                plateTypeSection = true;
            } 
            else if (qName.equals("Sublayout") && attributes.getValue("id") == null && sublayoutSection == false) 
            {
                sublayoutInfo = new SublayoutInfo();
                sublayoutSection = true;
            } 
            else if (qName.equals("Stack") && attributes.getValue("id") == null && stackPlanesSection == false) 
            {
                planeEntry = new StackPlaneEntry();
                stackPlanesSection = true;
            } 
            else if (qName.equals("Exposure") && attributes.getValue("id") == null && exposureSection == false) 
            {
                exposureInfo = new ExposureInfo();
                exposureSection = true;
            } 
            else if (qName.equals("Camera") && attributes.getValue("id") == null && cameraSection == false) 
            {
                cameraInfo = new CameraInfo();
                cameraSection = true;
            } 
            else if (qName.equals("SlowKinetic") && attributes.getValue("id") == null && slowKineticSection == false) 
            {
                slowKineticInfo = new SlowKineticInfo();
                slowKineticSection = true;
            } 
            else if (qName.equals("Objective") && attributes.getValue("id") == null && objectiveSection == false && exposureSection == false) 
            {
                objectiveInfo = new ObjectiveInfo();
                objectiveSection = true;
            } 
            else if (qName.equals("Filterchanger") && attributes.getValue("id") == null && filterChangerSection == false) 
            {
                filterChangerInfo = new FilterChangerInfo();
                filterChangerSection = true;
            } 
            else if (qName.equals("Filter") && attributes.getValue("id") == null && filterSection == false) 
            {
                filterInfo = new FilterInfo();
                filterSection = true;
            } 
            else if (qName.equals("DiscreteWavelength") && attributes.getValue("id") == null && lightSourceSection == false) 
            {
                lightSourceInfo = new LightSourceInfo();
                lightSourceSection = true;
            } 
            else if (qName.equals("Output") && attributes.getValue("id") == null && lightSourceOutputSection == false && lightSourceSection == true) 
            {
                lightSourceOutputInfo = new LightSourceOutputInfo();
                lightSourceOutputSection = true;
            } 
            else if (qName.equals("Field") && attributes.getValue("id") == null && sublayoutSection == true) 
            {
                fieldEntry = new FieldEntry();
            } 
            else if (qName.equals("Record") && attributes.getValue("id") == null && exposureSection == true) 
            {
                exposureRecordEntry = new ExposureRecordEntry();
            } 
            else if (qName.equals("Repetition") && attributes.getValue("id") == null && slowKineticSection == true) 
            {
                repetition = new Repetition();
            } 
            else if (qName.equals("Well") && attributes.getValue("id") == null) 
            {
                wellEntry = new WellEntry();
            } 
        }

        @Override
        public void endElement(String uri, String localName, String qName) 
        {
            String value = currentValue.toString();
            
            if ("Serial".equals(currentName)) {
                serial = value;
            }
            else if ("MeasurementID".equals(currentName)) {
                measurementID = value;
            }
            else if ("CurrentTemperature".equals(currentName)) {
                final double temp = Double.parseDouble(value);
                currentTemperature = new Temperature(temp, UNITS.CELSIUS);
            }
            else if ("CurrentCO2".equals(currentName)) {
                final Float co2;
                if (Float.parseFloat(value) > new Float(1.0)) {
                    co2 = Float.parseFloat(value) / 100;
                } else {
                    co2 = Float.parseFloat(value);
                }
                currentCO2 = new PercentFraction(co2);
            }
            else if ((objectiveSection == false) && ("Magnification".equals(currentName))) {
                magnification = Double.parseDouble(value);
            }
            else if (imageInfoSection != false) {
                if ("Row".equals(currentName)) {
                    imageInfo.row = Integer.parseInt(value);
                }
                else if ("Col".equals(currentName)) {
                    imageInfo.column = Integer.parseInt(value);
                }
                else if ("FieldID".equals(currentName)) {
                    imageInfo.field = Integer.parseInt(value);
                }
                else if ("PlaneID".equals(currentName)) {
                    imageInfo.plane = Integer.parseInt(value);
                }
                else if ("SlowKineticID".equals(currentName)) {
                    imageInfo.kinetic = Integer.parseInt(value);
                }
                else if ("ExposureID".equals(currentName)) {
                    imageInfo.exposure = Integer.parseInt(value);
                }
                else if ("RecordID".equals(currentName)) {
                    imageInfo.record = Integer.parseInt(value);
                }
                else if ("ChannelID".equals(currentName)) {
                    imageInfo.channel = Integer.parseInt(value);
                }
                else if ("FlimID".equals(currentName)) {
                    imageInfo.flim = Integer.parseInt(value);
                }
                else if ("AbsTime".equals(currentName)) {
                    imageInfo.absTime = value;
                }
            }
            else if (sublayoutSection != false) {
                if (fieldEntry != null) {
                    if ("X".equals(currentName)) {
                        fieldEntry.x = Double.parseDouble(value);
                    } 
                    else if ("Y".equals(currentName)) {
                        fieldEntry.y = Double.parseDouble(value);
                    }
                } 
                else if ("Plate".equals(currentName)) {
                    sublayoutInfo.plateName = value;
                } 
                else if ("Lens".equals(currentName)) {
                    sublayoutInfo.lensName = value;
                } 
                else if ("OverlapX".equals(currentName)) {
                    sublayoutInfo.overlapX = Double.parseDouble(value);
                } 
                else if ("OverlapY".equals(currentName)) {
                    sublayoutInfo.overlapY = Double.parseDouble(value);
                }
            }
            else if (plateTypeSection != false) {
                if ("ID".equals(currentName)) {
                    plateInfo.plateId = value;
                }
                else if ("Name".equals(currentName)) {
                    plateInfo.plateDescription = value;
                }
                else if ("FootPrintWidthX".equals(currentName)) {
                    plateInfo.footprint_width_x = Double.parseDouble(value);
                }
                else if ("FootPrintWidthY".equals(currentName)) {
                    plateInfo.footprint_width_y = Double.parseDouble(value);
                }
                else if ("StartX".equals(currentName)) {
                    plateInfo.start_x = Double.parseDouble(value);
                }
                else if ("StartY".equals(currentName)) {
                    plateInfo.start_y = Double.parseDouble(value);
                }
                else if ("EndX".equals(currentName)) {
                    plateInfo.end_x = Double.parseDouble(value);
                }
                else if ("EndY".equals(currentName)) {
                    plateInfo.end_y = Double.parseDouble(value);
                }
                else if ("StartRow".equals(currentName)) {
                    plateInfo.start_row = Integer.parseInt(value);
                }
                else if ("StartCol".equals(currentName)) {
                    plateInfo.start_column = Integer.parseInt(value);
                }
                else if ("EndRow".equals(currentName)) {
                    plateInfo.end_row = Integer.parseInt(value);
                }
                else if ("EndCol".equals(currentName)) {
                    plateInfo.end_column = Integer.parseInt(value);
                }
                else if ("WellFormTop".equals(currentName)) {
                    plateInfo.wellFormTop = value;
                }
                else if ("WellFormBottom".equals(currentName)) {
                    plateInfo.wellFormBottom = value;
                }
            }
            else if (wellEntry != null) {
                if ("Row".equals(currentName)) {
                    wellEntry.row = Integer.parseInt(value);
                }
                else if ("Col".equals(currentName)) {
                    wellEntry.column = Integer.parseInt(value);
                }
            }
            else if (exposureSection != false) {
                if ("Objective".equals(currentName)) {
                    exposureInfo.objectiveName = value;
                } else if ("Confocal".equals(currentName)) {
                    exposureInfo.confocal = Boolean.parseBoolean(value);
                } 
                if (exposureRecordEntry != null) {
                    if ("Name".equals(currentName)) {
                        exposureRecordEntry.channelName = value;
                    } 
                    else if ("FilterID".equals(currentName)) {
                        exposureRecordEntry.emissionFilterName = value;
                    } 
                    else if ("FilterName".equals(currentName)) {
                        exposureRecordEntry.emissionFilterWavelengths = value;
                    } 
                    else if ("BinningX".equals(currentName)) {
                        exposureRecordEntry.binningx = Integer.parseInt(value);
                    } 
                    else if ("BinningY".equals(currentName)) {
                        exposureRecordEntry.binningy = Integer.parseInt(value);
                    } 
                    else if ("Exposuretime".equals(currentName)) {
                        exposureRecordEntry.exposureTime = Double.parseDouble(value);
                    } 
                    else if ("Power".equals(currentName)) {
                        exposureRecordEntry.lightSourceOutputPower = Integer.parseInt(value);
                    } 
                    else if ("OutputID".equals(currentName)) {
                        exposureRecordEntry.lightSourceOutputName = value;
                    } 
                    else if ("OutputName".equals(currentName)) {
                        exposureRecordEntry.lightSourceWavelengths = value;
                    } 
                    else if ("Height".equals(currentName)) {
                        exposureRecordEntry.zHeight = Double.parseDouble(value);
                    } 
                    else if ("FilterchangerID".equals(currentName)) {
                        exposureRecordEntry.emissionFilterID = value;
                    } 
                    else if ("LightsourceID".equals(currentName)) {
                        exposureRecordEntry.lightSourceID = value;
                    } 
                }
            } 
            else if (slowKineticSection != false) {
                if ("Combined".equals(currentName)) {
                    slowKineticInfo.combined = Boolean.parseBoolean(value);
                } else if ("Alarm".equals(currentName)) {
                    slowKineticInfo.alarm = Boolean.parseBoolean(value);
                } else if ("UseBreak".equals(currentName)) { 
                    slowKineticInfo.usebreak = Boolean.parseBoolean(value);
                } 
                if (repetition != null) {
                    if ("TimepointNr".equals(currentName)) {
                        repetition.timepoint = Integer.parseInt(value);
                    } else if ("StartDelay".equals(currentName)) {
                        repetition.startdelay = value;
                    } else if ("Pause".equals(currentName)) {
                        repetition.pause = Boolean.parseBoolean(value);
                    } else if ("Fast".equals(currentName)) {
                        repetition.fast = Boolean.parseBoolean(value);
                    }
                }
            }
            else if (cameraSection != false) {
                if ("PixelsX".equals(currentName)) {
                    cameraInfo.xpixels = Integer.parseInt(value);
                } else if ("PixelSizeX".equals(currentName)) {
                    cameraInfo.xpixelsize = Double.parseDouble(value);
                } else if ("PixelsY".equals(currentName)) { 
                    cameraInfo.ypixels = Integer.parseInt(value);
                } else if ("PixelSizeY".equals(currentName)) { 
                    cameraInfo.ypixelsize = Double.parseDouble(value);
                } 
            }
            else if (objectiveSection != false) {
                if ("ID".equals(currentName)) {
                    objectiveInfo.ID = value;
                } else if ("Name".equals(currentName)) {
                    objectiveInfo.name = value;
                } else if ("Magnification".equals(currentName)) { 
                    objectiveInfo.magnification = Integer.parseInt(value);
                }
            }
            else if (filterChangerSection != false && filterSection == false) {
                if ("ID".equals(currentName)) {
                    filterChangerInfo.ID = value;
                } else if ("Name".equals(currentName)) {
                    filterChangerInfo.name = value;
                } 
            }
            else if (filterSection != false) {
                if ("ID".equals(currentName)) {
                    filterInfo.ID = value;
                } else if ("Name".equals(currentName)) {
                    filterInfo.range = value;
                } else if ("MainWavelength".equals(currentName)) {
                    filterInfo.wavelength = Integer.parseInt(value);
                } 
            }
            else if (lightSourceSection != false && lightSourceOutputSection == false) {
                if ("ID".equals(currentName)) {
                    lightSourceInfo.ID = value;
                } else if ("Name".equals(currentName)) {
                    lightSourceInfo.type = value;
                }
            }
            else if (lightSourceOutputSection != false) {
                if ("ID".equals(currentName)) {
                    lightSourceOutputInfo.ID = value;
                } else if ("Name".equals(currentName)) {
                    lightSourceOutputInfo.range = value;
                } else if ("MainWavelength".equals(currentName)) {
                    lightSourceOutputInfo.wavelength = Integer.parseInt(value);
                } 
            }
            
            currentName = null;

            if (qName.equals("Field") && sublayoutInfo != null && sublayoutSection == true) {
                sublayoutInfo.fieldEntries.add(fieldEntry);
                fieldEntry = null;
            }
            else if (qName.equals("Well") && wellEntry != null) {
                wells.add(wellEntry);
                wellEntry = null;
            }
            else if (qName.equals("Stack") && planeEntry != null) {
                planes.add(planeEntry);
                planeEntry = null;
            }
            else if (qName.equals("Repetition") && repetition != null) {
                slowKineticInfo.timepoints.add(repetition);
                repetition = null;
            }
            else if (qName.equals("Record") && exposureRecordEntry != null && exposureSection == true) {
                exposureInfo.exposureRecordEntries.add(exposureRecordEntry);
                exposureRecordEntry = null;
            } 
            else if (qName.equals("Objective") && objectiveInfo != null && objectiveSection == true) {
                objectives.add(objectiveInfo);
                objectiveInfo = null;
                objectiveSection = false;
            }
            else if (qName.equals("Filter") && filterInfo != null && filterSection == true) {
                filterChangerInfo.filters.add(filterInfo);
                filterInfo = null;
                filterSection = false;
            }
            else if (qName.equals("Output") && lightSourceOutputInfo != null && lightSourceOutputSection == true) {
                lightSourceInfo.outputs.add(lightSourceOutputInfo);
                lightSourceOutputInfo = null;
                lightSourceOutputSection = false;
            }
            else if (qName.equals("Filterchanger") && filterChangerInfo != null && filterChangerSection == true) {
                filterChangers.add(filterChangerInfo);
                filterChangerInfo = null;
                filterChangerSection = false;
            }
            else if (qName.equals("DiscreteWavelength") && lightSourceInfo != null && lightSourceSection == true) {
                lightSources.add(lightSourceInfo);
                lightSourceInfo = null;
                lightSourceSection = false;
            }
            else if (qName.equals("PlateType") && plateInfo != null && plateTypeSection == true) {
                plateTypeSection = false;
            }
            else if (qName.equals("Sublayout") && sublayoutInfo != null && sublayoutSection == true) {
                sublayoutSection = false;
            }
            else if (qName.equals("Image") && imageInfo != null && imageInfoSection == true) {
                imageInfoSection = false;
            }
            else if (qName.equals("Exposure") && exposureInfo != null && exposureSection == true) {
                exposureSection = false;
            }
            else if (qName.equals("Camera") && cameraInfo != null && cameraSection == true) {
                cameraSection = false;
            }
            else if (qName.equals("SlowKinetic") && slowKineticInfo != null && slowKineticSection == true) {
                slowKineticSection = false;
            }
        }
    }

    class ImageInfo {
        public Integer row;
        public Integer column;
        public Integer field;
        public Integer plane;
        public Integer kinetic;
        public Integer exposure;
        public Integer record;
        public Integer channel;
        public Integer flim;
        public String absTime;
    }
    
    class PlateInfo {
        public String plateId;      
        public String plateDescription;
        public String wellFormTop;
        public String wellFormBottom;

        public int start_column;
        public int end_column;
        public int start_row;
        public int end_row;
        
        public double start_x;
        public double end_x;
        public double start_y;
        public double end_y;
        public double footprint_width_x;
        public double footprint_width_y;
    }

    class SublayoutInfo {
        public ArrayList<FieldEntry> fieldEntries = new ArrayList<FieldEntry>();
        public String lensName;
        public String plateName;
        public double overlapX;
        public double overlapY;
    }

    class FieldEntry {
        public double x;
        public double y;
    }

    class StackPlaneEntry {
        public double z;
        public int exposureId;
    }

    class WellEntry {
        public Integer column;
        public Integer row;
    }
    
    class ExposureInfo {
        public ArrayList<ExposureRecordEntry> exposureRecordEntries = new ArrayList<ExposureRecordEntry>();
        public String objectiveName;
        public Boolean confocal;
    }
    
    class ExposureRecordEntry {
        public Double zHeight;
        public String emissionFilterID;
        public String emissionFilterName;
        public String emissionFilterWavelengths;
        public String channelName;
        public Double exposureTime;
        public String lightSourceID;
        public String lightSourceOutputName;
        public String lightSourceWavelengths;
        public Integer lightSourceOutputPower;
        public Integer binningx;
        public Integer binningy;
    }
    
    class CameraInfo {
        public Integer xpixels;
        public Integer ypixels;
        public Double xpixelsize;
        public Double ypixelsize;
    }
    
    class SlowKineticInfo {
        public ArrayList<Repetition> timepoints = new ArrayList<Repetition>();
        public Boolean combined;
        public Boolean alarm;
        public Boolean usebreak;
    }
    
    class Repetition {
        public Integer timepoint;
        public String startdelay;
        public Boolean pause;
        public Boolean fast;
    }
    
    class ObjectiveInfo {
        public String ID;
        public String name;
        public Integer magnification;
    }
    
    class LightSourceInfo {
        public String ID;
        public String type;
        public ArrayList<LightSourceOutputInfo> outputs = new ArrayList<LightSourceOutputInfo>();
    }
    
    class LightSourceOutputInfo {
        public String ID;
        public String range;
        public Integer wavelength;
    }
    
    class FilterChangerInfo {
        public String ID;
        public String name;
        public ArrayList<FilterInfo> filters = new ArrayList<FilterInfo>();
    }
    
    class FilterInfo {
        public String ID;
        public String name;
        public String range;
        public Integer wavelength;
    }
}
