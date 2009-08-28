/* Mesquite source code.  Copyright 1997-2009 W. Maddison and D. Maddison. Version 2.7, August 2009.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html) *//* * XMLConstants.java * * Created on June 20, 2003, 8:09 AM */package mesquite.tol.lib;/** * Constants interface that has the names of XML tags and attributes for the ToL Project * @author dmandel * */public interface XMLConstants {    public static final String OTHER_ACTIVE_DOWNLOAD = "OTHER_ACTIVE_DOWNLOAD";    public static final String CHARSET_NAME = "ISO-8859-1";    public static final String DONTPUBLISH = "DONTPUBLISH";    public static final String DONTPUBLISHCHANGED = "DONTPUBLISHCHANGED";    	public static final String ISEDITOR = "ISEDITOR";    public static final String DOWNLOADS = "DOWNLOADS";    public static final String DOWNLOAD_ID = "DOWNLOAD_ID";    public static final String TREESTRUCTURE = "TREESTRUCTURE";    public static final String FILENAME = "FILENAME";    public static final String NODEID = "NODEID";    public static final String NODES = "NODES";    public static final String NODE = "NODE";    public static final String NAME = "NAME";    public static final String NAMECOMMENT = "NAMECOMMENT";    public static final String FIRSTNAME = "FIRSTNAME";    public static final String LASTNAME = "LASTNAME";    public static final String NAMECHANGED = "NAMECHANGED";    public static final String TRUE = "TRUE";    public static final String FALSE = "FALSE";    public static final String ID = "ID";    public static final String LOCKED = "LOCKED";    public static final String ONE = "1";    public static final String ZERO = "0";    public static final String LOCK_INFO = "LOCK_INFO";    public static final String TIMESTAMP = "TIMESTAMP";    public static final String IPADDR = "IPADDR";    public static final String DATE_TIME = "DATE_TIME";    public static final String USER = "USER";    public static final String TYPE = "TYPE";    public static final String EXTINCT = "EXTINCT";    public static final String EXTINCTCHANGED = "EXTINCTCHANGED";    public static final String CONFIDENCE = "CONFIDENCE";    public static final String CONFIDENCECHANGED = "CONFIDENCECHANGED";    public static final String PHYLESIS = "PHYLESIS";    public static final String PHYLESISCHANGED = "PHYLESISCHANGED";    public static final String LEAF = "LEAF";    public static final String LEAFCHANGED = "LEAFCHANGED";    public static final String NODERANK = "NODERANK";    public static final String NODERANKCHANGED = "NODERANKCHANGED";    public static final String PAGE = "PAGE";    public static final String PAGES = "PAGES";    public static final String PAGECHANGED = "PAGECHANGED";    public static final String DESCRIPTION = "DESCRIPTION";    public static final String DESCRIPTIONCHANGED = "DESCRIPTIONCHANGED";    public static final String OTHERNAMES = "OTHERNAMES";    public static final String OTHERNAMESCHANGED = "OTHERNAMESCHANGED";    public static final String OTHERNAME = "OTHERNAME";    public static final String ISLABEL = "ISLABEL";    public static final String CHILDRENCHANGED = "CHILDRENCHANGED";    public static final String FIRSTONLINECHANGED = "FIRSTONLINECHANGED";    public static final String FIRSTONLINE = "FIRSTONLINE";    public static final String OPTIONSCHANGED = "OPTIONSCHANGED";    public static final String OPTIONS = "OPTIONS";    public static final String STATUS = "STATUS";    public static final String STATUSCHANGED = "STATUSCHANGED";    public static final String WRITEASLIST = "WRITEASLIST";    public static final String WRITECHANGED = "WRITECHANGED";    public static final String GENBANKCHANGED = "GENBANKCHANGED";    public static final String GENBANK = "GENBANK";    public static final String TREEBASECHANGED = "TREEBASECHANGED";    public static final String TREEBASE = "TREEBASE";    public static final String TITLECHANGED = "TITLECHANGED";    public static final String TITLE = "TITLE";    public static final String SUBTITLECHANGED = "SUBTITLECHANGED";    public static final String SUBTITLE = "SUBTITLE";    public static final String PRETREETEXTCHANGED = "PRETREETEXTCHANGED";    public static final String PRETREETEXT = "PRETREETEXT";    public static final String POSTTREETEXTCHANGED = "POSTTREETEXTCHANGED";    public static final String POSTTREETEXT = "POSTTREETEXT";    public static final String IMGCAPTIONCHANGED = "IMGCAPTIONCHANGED";    public static final String IMGCAPTION = "IMGCAPTION";    public static final String TAGGED_IMGCAPTIONCHANGED = "TAGGED_IMGCAPTIONCHANGED";    public static final String TAGGED_IMGCAPTION = "TAGGED_IMGCAPTION";    public static final String ACKSCHANGED = "ACKSCHANGED";    public static final String ACKS = "ACKS";    public static final String CONTENTCHANGEDDATECHANGED = "CONTENTCHANGEDDATECHANGED";    public static final String CONTENTCHANGEDDATE = "CONTENTCHANGEDDATE";    public static final String COPYRIGHTCHANGED = "COPYRIGHTCHANGED";    public static final String COPYRIGHT = "COPYRIGHT";    public static final String COPYRIGHTOWNER = "COPYRIGHTOWNER";    public static final String COPYRIGHTURL = "COPYRIGHTURL";    public static final String COPYRIGHTEMAIL = "COPYRIGHTEMAIL";    public static final String COPYRIGHTDATE = "COPYRIGHTDATE";        public static final String DATE = "DATE";    public static final String DATECHANGED = "DATECHANGED";    public static final String HOLDER = "HOLDER";    public static final String HOLDERCHANGED = "HOLDERCHANGED";    public static final String IMAGECHANGED = "IMAGECHANGED";    public static final String IMAGELIST = "IMAGELIST";    public static final String IMAGE = "IMAGE";    public static final String IMAGES = "IMAGES";       public static final String FOSSIL = "FOSSIL";    public static final String SEX = "SEX";    public static final String COLLECTIONDATE = "COLLECTIONDATE";    public static final String COLLECTION = "COLLECTION";    public static final String COLLECTIONAC = "COLLECTIONAC";    public static final String LOCATION = "LOCATION";    public static final String GEOLOCATION = "GEOLOCATION";    public static final String CAPTION = "CAPTION";    public static final String SCIENTIFICNAME = "SCIENTIFICNAME";    public static final String SCIENTIFICNAMETWO = "SCIENTIFICNAMETWO";    public static final String COMMONNAME = "COMMONNAME";    public static final String COMMONNAMETWO = "COMMONNAMETWO";        public static final String CREATOR = "CREATOR";    public static final String IDENTIFIER = "IDENTIFIER";    public static final String STAGE = "STAGE";    public static final String BODYPART = "BODYPART";    public static final String SIZE = "SIZE";    public static final String VIEW = "VIEW";    public static final String PERIOD = "PERIOD";    public static final String COLLECTOR = "COLLECTOR";        public static final String URL = "URL";    public static final String USECONTENT = "USECONTENT";    public static final String LINKED_MENU = "LINKED_MENU";    public static final String LINKED_NODE_ID = "LINKED_NODE_ID";    public static final String LINKED_ACC_ID = "LINKED_ACC_ID";    public static final String ORDER = "ORDER";    public static final String LISTCHANGED = "LISTCHANGED";    public static final String TEXTLIST = "TEXTLIST";    public static final String CHANGED = "CHANGED";    public static final String HEADING = "HEADING";    public static final String TEXT = "TEXT";   // public static final String SEQUENCE = "SEQUENCE";    public static final String ACCESSORYCHANGED = "ACCESSORYCHANGED";    public static final String ACCESSORYPAGES = "ACCESSORYPAGES";    public static final String ACCESSORYPAGE = "ACCESSORYPAGE";    public static final String LINKCHANGED = "LINKCHANGED";    public static final String LINK = "LINK";        public static final String LINKS = "LINKS";    public static final String AUTHORCHANGED = "AUTHORCHANGED";    public static final String AUTHORLIST = "AUTHORLIST";    public static final String FULLNAME = "FULLNAME";    public static final String EMAIL = "EMAIL";     public static final String CONTACT = "CONTACT";    public static final String AUTHOR = "AUTHOR";    public static final String IS_AUTHOR = "IS_AUTHOR";    public static final String ADDRESS = "ADDRESS";    public static final String ERROR = "ERROR";    public static final String ERRORNUM = "ERRORNUM";    public static final String ERRORTEXT = "ERRORTEXT";    public static final String ANCESTORS = "ANCESTORS";    public static final String ANCESTORS_INFO = "ANCESTORS_INFO";    public static final String ANCESTOR_INFO = "ANCESTOR_INFO";    public static final String TREE = "TREE";    public static final String DEPTH = "DEPTH";    public static final String HASPAGE = "HASPAGE";    public static final String ANCESTORPAGEID = "ANCESTORWITHPAGE";    public static final String FILECHANGED = "FILECHANGED";    public static final String UNLOCKED = "UNLOCKED";       public static final String SEQUENCECHANGED = "SEQUENCECHANGED";    public static final String CHILDCOUNT = "CHILDCOUNT";    public static final String IMAGEURL = "IMAGEURL";    public static final String TEXTSECTION = "TEXTSECTION";    public static final String INDENT = "INDENT";    public static final String REFERENCES = "REFERENCES";    public static final String REFERENCESCHANGED = "REFERENCESCHANGED";    public static final String INTERNETLINKS = "INTERNETLINKS";    public static final String INTERNETLINKSCHANGED = "INTERNETLINKSCHANGED";       public static final String USERNAME = "USERNAME";    public static final String SETTINGS = "SETTINGS";    public static final String CUSTOM_CURSORS = "CUSTOM_CURSORS";    public static final String EDITORJAR_TIMESTAMP = "EDITORJAR_TIMESTAMP";    public static final String SUPPORTJAR_TIMESTAMP = "SUPPORTJAR_TIMESTAMP";       public static final String MESSAGE = "message";    public static final String MESSAGENAME = "name";    public static final String MESSAGETEXT = "text";    public static final String SKELETAL = "Skeletal";    public static final String TEMPORARY = "Temporary";    public static final String UNDER_CONSTRUCTION = "Under Construction";    public static final String COMPLETE = "Complete";    public static final String TOL_REVIEWED = "ToL Reviewed";    public static final String PEER_REVIEWED = "Peer Reviewed";    public static final String MENU = "MENU";    public static final String MENUONSERVER = "MENUONSERVER";    public static final String PAGETITLE = "PAGETITLE";    public static final String AUTHORS = "AUTHORS";    public static final String AUTHDATE = "AUTHDATE";    public static final String SHOW = "SHOW";    public static final String SHOWAUTHORITYCHANGED = "SHOWAUTHORITYCHANGED";    public static final String SHOWAUTHORITY = "SHOWAUTHORITY";    public static final String SHOWAUTHORITYCONTAINING = "SHOWAUTHORITYCONTAINING";    public static final String SHOWPREFAUTHORITYCHANGED = "SHOWPREFAUTHORITYCHANGED";    public static final String SHOWPREFAUTHORITY = "SHOWPREFAUTHORITY";    public static final String SHOWIMPAUTHORITYCHANGED = "SHOWIMPAUTHORITYCHANGED";    public static final String SHOWIMPAUTHORITY = "SHOWIMPAUTHORITY";    public static final String AUTHORITY = "AUTHORITY";    public static final String ISPREFERRED = "ISPREFERRED";    public static final String ISIMPORTANT = "ISIMPORTANT";    public static final String MATCH = "MATCH";    public static final String MATCHES = "MATCHES";    public static final String ROOT = "ROOT";    public static final String MODIFIEDDATE = "MODIFIEDDATE";    public static final String UPLOADDATE = "UPLOADDATE";    public static final String DOWNLOADDATE = "DOWNLOADDATE";    public static final String CHECKED_OUT_FILE = "CHECKED_OUT_FILE";    public static final String UPLOAD_BATCH = "UPLOAD_BATCH";    public static final String UPLOAD_BATCHES = "UPLOAD_BATCHES";        public static final String BATCHID = "BATCHID";    public static final String SUBMITTED = "SUBMITTED";    public static final String DOWNLOADED = "DOWNLOADED";    public static final String LAST_USER = "LAST_USER";    public static final String LAST_DATE = "LAST_DATE";    public static final String ROOT_GROUP = "ROOT_GROUP";        public static final String ROOT_GROUP_ID = "ROOT_GROUP_ID";        public static final String ROOTNODE_ID = "ROOTNODE_ID";    public static final String PERMISSIONS = "PERMISSIONS";    public static final String PERMISSION = "PERMISSION";        public static final String SUCCESS = "SUCCESS";    public static final String FAILURE = "FAILURE";        public static final String UNDOABLE_UPLOAD = "UNDOABLE_UPLOAD";    public static final String usermessages = "usermessages";    public static final String NOID = "NOID";    public static final String WRONG_PASSWORD = "WRONG_PASSWORD";    public static final String PASSWORD = "PASSWORD";    public static final String CAN_PUSH_PUBLIC = "CAN_PUSH_PUBLIC";    public static final String IS_SOLE_AUTHOR = "IS_SOLE_AUTHOR";    public static final String NEW_FROM_SERVER = "NEW_FROM_SERVER";    public static final String PARENTPAGE_NAME = "PARENTPAGE_NAME";    public static final String PARENTPAGE_NODE_ID = "PARENTPAGE_NODE_ID";    public static final String CONTRIBUTOR = "CONTRIBUTOR";    public static final String ORGANIZATION = "ORGANIZATION";    public static final String HOMEPAGE = "HOMEPAGE";    public static final String CONTRIBUTORLIST = "CONTRIBUTORLIST";    public static final String CONTRIBUTORCHANGED = "CONTRIBUTORCHANGED";    public static final String WIDTH = "WIDTH";    public static final String HEIGHT = "HEIGHT";    public static final String IMAGEID = "IMAGEID";    public static final String REMOVED = "REMOVED";    public static final String ALTTEXT = "ALTTEXT";    public static final String COMMENTS = "COMMENTS";    public static final String DONTSHOWCLOSE = "DONTSHOWCLOSE";    public static final String PUBLICDOMAIN = "PUBLICDOMAIN";    public static final String PREFERENCES = "PREFERENCES";    public static final String PREFERENCE = "PREFERENCE";       public static final String PRINTIMAGEDATA = "PRINTIMAGEDATA";    public static final String PRINTIMAGECAPTION = "PRINTIMAGECAPTION";    public static final String INSTITUTION = "INSTITUTION";    public static final String NOTES = "NOTES";    public static final String NOEMAIL = "NOEMAIL";    public static final String NOADDRESS = "NOADDRESS";    public static final String PHONE = "PHONE";    public static final String FAX = "FAX";    public static final String SPECIMEN = "SPECIMEN";    public static final String BODYPARTS = "BODYPARTS";    public static final String ULTRASTRUCTURE = "ULTRASTRUCTURE";    public static final String HABITAT = "HABITAT";    public static final String EQUIPMENT = "EQUIPMENT";    public static final String PEOPLE = "PEOPLE";    public static final String CREATIONDATE = "CREATIONDATE";    public static final String ALIVE = "ALIVE";    public static final String BEHAVIOR = "BEHAVIOR";    public static final String VOUCHERNUMBER = "VOUCHERNUMBER";    public static final String VOUCHERNUMBERCOLLECTION = "VOUCHERNUMBERCOLLECTION";    public static final String IMAGETYPE = "IMAGETYPE";    public static final String UPLOADRESULTS = "UPLOADRESULTS";    public static final String VERSIONS = "VERSIONS";    public static final String VERSION = "VERSION";    public static final String IS_ARTICLE = "IS_ARTICLE";    public static final String IS_TREEHOUSE = "IS_TREEHOUSE";    public static final String NEW_VERSION = "NEW_VERSION";    public static final String PAGEREMOVED = "PAGEREMOVED";    public static final String PAGEADDED = "PAGEADDED";    public static final String INCOMPLETESUBGROUPS = "INCOMPLETESUBGROUPS";    public static final String INCOMPLETESUBGROUPSCHANGED = "INCOMPLETESUBGROUPSCHANGED";    public static final String PRIORITY = "PRIORITY";    public static final String PRIORITYCHANGED = "PRIORITYCHANGED";    public static final String CHARACTER = "CHARACTER";    public static final String TAXON = "TAXON";    public static final String NEWID = "NEWID";    public static final String image = "image";    public static final String filename = "filename";    public static final String copyrightowner = "copyrightowner";    public static final String copyrightemail = "copyrightemail";    public static final String copyrighturl = "copyrighturl";    public static final String copyrightdate = "copyrightdate";    public static final String license = "license";    public static final String restricted = "1";    public static final String tolusenomod = "2n";    public static final String tolusemod = "2y";    public static final String tolsharenomod = "3n";    public static final String tolsharemod = "3y";    public static final String cc = "4";    public static final String pd = "pd";    public static final String reference = "reference";    public static final String creator = "creator";    public static final String acknowledgements = "acknowledgements";    public static final String group = "group";    public static final String y = "y";    public static final String n = "n";    public static final String specimen = "specimen";    public static final String bodyparts = "bodyparts";    public static final String ultrastructure = "ultrastructure";    public static final String habitat = "habitat";    public static final String equipment = "equipment";    public static final String people = "people";    public static final String subject = "subject";    public static final String keywords = "keywords";    public static final String geolocation = "geolocation";    public static final String time = "time";    public static final String condition = "condition";    public static final String live = "live";    public static final String l = "l";    public static final String dead = "dead";    public static final String d = "d";    public static final String model = "model";    public static final String m = "m";    public static final String fossil = "fossil";    public static final String f = "f";    public static final String period = "period";    public static final String scientificname = "scientificname";    public static final String identifier = "identifier";    public static final String behavior = "behavior";    public static final String sex = "sex";    public static final String stage = "stage";    public static final String partofbody = "partofbody";    public static final String view = "view";    public static final String size = "size";    public static final String collection = "collection";    public static final String type = "type";    public static final String vouchernumber = "vouchernumber";    public static final String vouchercollection = "vouchercollection";    public static final String collector = "collector";    public static final String comments = "comments";    public static final String alt = "alt";    public static final String imagetype = "imagetype";    public static final String artistic = "artistic";    public static final String technical = "technical";    public static final String notes = "notes";    public static final String photo = "photo";    public static final String painting = "painting";    public static final String drawing = "drawing";    public static final String diagram = "diagram";    public static final String person = "person";    public static final String surname = "surname";    public static final String firstname = "firstname";    public static final String email = "email";    public static final String dontpublishemail = "dontpublishemail";    public static final String webpageurl = "webpageurl";    public static final String institution = "institution";    public static final String address = "address";    public static final String bio = "bio";    public static final String grouphaspermission = "grouphaspermission";    public static final String groupnopermission = "groupnopermission";    public static final String coordinator = "coordinator";    public static final String phone = "phone";    public static final String fax = "fax";	public static final String interests = "interests";	public static final String geographicareainterest = "geographicareainterest";	public static final String otherinterests = "otherinterests";	public static final String permission = "permission"; 	public static final String sendemail = "sendemail";	public static final String direction = "direction";	public static final String genename = "genename";	public static final String sequence = "sequence";	public static final String branch = "branch";	public static final String attribute = "attribute";	public static final String name = "name";	public static final String value = "value";	public static final String PARENTGROUP = "PARENTGROUP";	public static final String COUNT = "COUNT";	public static final String OBSOLETEMESSAGE = "OBSOLETEMESSAGE";	public static final String VMVERSION = "VMVERSION";	public static final String ITALICIZE_NAME = "ITALICIZENAME";	public static final String IS_NEW_COMBINATION = "IS_NEW_COMBINATION";	public static final String COMBINATION_AUTHOR = "COMBINATION_AUTHOR";	public static final String COMBINATION_DATE = "COMBINATION_DATE";	public static final String SOURCE_DB = "SOURCE_DB";	public static final String SOURCE_DB_ID = "SOURCE_DB_ID";	}