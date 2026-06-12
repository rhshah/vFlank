import csv, pysam, ConfigParser
#from primerClasses import sequence, FusionSequence

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Emboss.Applications import Primer3Commandline



def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

def get_fusion_partner(chr, pos, length, fasta, type, direction):
    FASTA = pysam.Fastafile(fasta)
    sequence = ""
    try:
        if direction == "+":
            if type == 1:
                start = int(pos) + int(length)
                end = pos
            elif type ==2:
                start = pos
                end = int(pos) + int(length)

        elif direction == "-":
            if type == 1:
                start = int(pos) - int(length)
                end = pos
            elif type ==2:
                end = int(pos) - int(length)
                start = pos
        print "obtaining sequence for %s, %s, %s, %s"%(chr, start, end, fasta)
        sequence = Seq(str(FASTA.fetch(region="%s:%s-%s" %(chr, start, end))))
                           
    except Exception as e:
        print e
    if direction == "+":
        return str(sequence)
    elif direction == "-":
        return str(sequence.reverse_complement())
    
if __name__ == '__main__':
       
    import sys
    Config = ConfigParser.ConfigParser()
    #Config.read("/dmp/hot/zehira/ctDNA_primers/config.cfg")
    if(sys.argv[1] == ""):
	print "Please provide the config file."
	sys.exit(1)
    if(sys.argv[2] == ""):
        print "Please provide the direction of transcript for 1st break point 1=+,2=-."
        sys.exit(1)
    if(sys.argv[3] == ""):
        print "Please provide the direction of transcript for 2nd break point 1=+,2=-."
        sys.exit(1)
    if(sys.argv[4] == ""):
        print "Please provide name of the output fasta file"
        sys.exit(1)
    Config.read(sys.argv[1])
    FusionOptions = ConfigSectionMap("FusionOptions")
    ProbeOptions = ConfigSectionMap("ProbeOptions")
    GeneralOptions= ConfigSectionMap("GeneralOptions")
    fusion_half_length = int(FusionOptions['fusion_length'])/2
    
    fusion_seq = SeqRecord(Seq(get_fusion_partner(FusionOptions['partner_1_chr'], FusionOptions['partner_1_start'], -fusion_half_length, GeneralOptions['fasta'], int(sys.argv[2]), FusionOptions['partner_1_dir']) ) + "-" + Seq(get_fusion_partner(FusionOptions['partner_2_chr'], FusionOptions['partner_2_start'], fusion_half_length, GeneralOptions['fasta'], int(sys.argv[3]), FusionOptions['partner_2_dir'])))
    #fusion_seq = SeqRecord(Seq("actagcatgcatgctagctagctagtcgatc"))
    fusion_seq.id = FusionOptions['patient_id']
    fusion_seq.name = FusionOptions['fusion_name'] + "_" + FusionOptions['partner_1_chr'] + ":" + FusionOptions['partner_1_start'] + "_" + FusionOptions['partner_2_chr'] + ":" + FusionOptions['partner_2_start']
    fusion_seq.description = FusionOptions['description']
    probe_feature = SeqFeature(FeatureLocation(190, 205), type="probe", strand=1)
    print "Fusion Sequence"
    print fusion_seq.format("fasta")
    
    print "\nFusion Sequnce Reverse compliment"
    a = fusion_seq.reverse_complement()
    print a.format("fasta")
    outfile = open(sys.argv[4], "w")
    SeqIO.write(fusion_seq, outfile, "fasta")
    outfile.close()
    try:
        probe = probe_feature.extract(fusion_seq)
        #sequence = get_fusion_partner(FusionOptions['partner_1_chr'], FusionOptions['partner_1_start'], -fusion_half_length, GeneralOptions['fasta'], 1) + "-" + get_fusion_partner(FusionOptions['partner_2_chr'], FusionOptions['partner_2_start'], fusion_half_length, GeneralOptions['fasta'], 2)
        primer_cl = Primer3Commandline(sequence = sys.argv[4], hybridprobe=True)
        primer_cl.osizeopt = ProbeOptions['probe_length']
        primer_cl.otmopt = 67
        primer_cl.ogcopt = 50
        primer_cl.outfile = "myresults.out"

        stdout, stderr = primer_cl()
        print "Command line for primer3: %s"%(primer_cl)
        print "STDout: %s\nSTDerr: %s"%(stdout, stderr)
    except Exception as e:
        print "Error: %s"%(e)
    
    
    
    
    


    
    




