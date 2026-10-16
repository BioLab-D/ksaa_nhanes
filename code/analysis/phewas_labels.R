# PheWAS 图表共用:健康类别的显示名与配色、变量的规范标签(SAS 标签多为缩写);未列出者用 tidy() 整理后的 SAS 标签
CLS_ORDER <- c("objective", "function", "symptom", "disease", "self_rated", "healthcare")
CLS_NAME <- c(objective = "Objective measures", function_ = "Physical function", symptom = "Symptoms", disease = "Disease history",
              self_rated = "Self-/physician-rated health", healthcare = "Health-care use"); names(CLS_NAME)[2] <- "function"
CLS_COL <- c(objective = "#2a78d6", `function` = "#8e5bd1", symptom = "#eb6834", disease = "#1baf7a", self_rated = "#eda100", healthcare = "#a3a29b")
tidy <- function(x) trimws(gsub("\\s+", " ", sub("\\s*\\((all ages|age ?[0-9]+\\+? ?(yrs|years)?|[0-9]+ (months|years|yrs) and over|[0-9]+\\+ ?(months|years)|age60\\+|age 60\\+|age5\\+|age 5\\+)\\)?\\s*$", "", x)))
NICE <- c(BDPTRBMD = "Trochanter BMD", BDPTOBMD = "Total hip BMD", BDPINBMD = "Intertrochanter BMD", BDPWTBMD = "Ward's triangle BMD", BDPFNBMD = "Femoral neck BMD",
  BDPTOARE = "Total hip bone area", PEP1 = "Abnormal locomotion (examination)", PEP1A1 = "Limp or shuffle (examination)", PEP1A2 = "Other gait abnormality (examination)",
  HAY11D = "Walks slowly or shuffles (observed)", HAY11C = "Uses cane, crutches or walker (observed)", HAY11G = "Speech problems (observed)",
  PFPNSTPA = "Steps to walk 8 ft", PFPNSTPB = "Steps to walk 8 ft (second trial)", PFPTWLKA = "8-ft walk time", PFPTWLKB = "8-ft walk time (second trial)",
  PFPRFLEX = "Limited right hip and knee flexion", PFPLFLEX = "Limited left hip and knee flexion", PFPDEVIC = "Walking aid used in the 8-ft walk", PFPTTAND = "Unable to hold tandem stand for 10 s",
  PFPUSTND = "Needs arms or unable to stand from chair", PFPTSTND = "Five chair stands time", PFPRSER = "Limited right shoulder external rotation", PFPLSIR = "Limited left shoulder internal rotation", PFPPAIN = "Pain on walking (examination)",
  HAH1 = "Difficulty walking a quarter mile", HAH2 = "Difficulty walking 10 steps without resting", HAH3 = "Difficulty stooping, crouching or kneeling", HAH4 = "Difficulty lifting or carrying 10 lb",
  HAH5 = "Difficulty doing household chores", HAH6 = "Difficulty preparing own meals", HAH7 = "Difficulty managing money", HAH8 = "Difficulty walking from room to room", HAH9 = "Difficulty standing up from an armless chair",
  HAH10 = "Difficulty getting in or out of bed", HAH11 = "Difficulty eating", HAH12 = "Difficulty dressing", HAH13 = "Needs help with personal care", HAH14 = "Needs help with routine needs", HAH15 = "Uses a cane or other walking aid",
  HAS24 = "Housework limited by a health problem", HAS25 = "Activities limited by a health problem",
  PEP13E1A = "Physician-estimated difficulty walking a quarter mile", PEP13E2A = "Physician-estimated difficulty running 100 yards", PEP13E3A = "Physician-estimated difficulty stooping",
  PEP13E4A = "Physician-estimated difficulty with small hand movements", PEP13E5A = "Physician-estimated difficulty with heavy housework", PEP13A = "Physician-rated health status", HAB1 = "Self-rated general health",
  HAG34 = "Chronic limp (more than 3 months)", HAG22 = "Knee pain on most days for at least 6 weeks", HAG25 = "Knee swelling on most days for at least 6 weeks", HAG27 = "Morning knee stiffness on most days for at least 6 weeks",
  HAG20 = "Morning hand stiffness on most days for at least 6 weeks", HAG14 = "At least one fall in the past 12 months", HAF14 = "Leg pain on walking", HAP17B = "Trials needed to learn three words",
  HAC1A = "Physician-diagnosed arthritis", HAG5A = "Self-reported hip fracture", HAG11 = "Physician-diagnosed osteoporosis", HAJ9 = "Physician-diagnosed gallstones", HAD1 = "Physician-diagnosed diabetes",
  PHPINSU = "Current insulin use", HAC5A1 = "Diabetes in grandmother", HAE2 = "Physician-diagnosed hypertension", MAPF25 = "Hysterectomy",
  APPSI = "Serum alkaline phosphatase", PEP10B4 = "Knee pain on passive motion (examination)", BMPSUB = "Subscapular skinfold", BMPTRI = "Triceps skinfold", BMPSUP = "Suprailiac skinfold", BMPTHI = "Thigh skinfold",
  BMPTHICI = "Thigh circumference", BMPARMC = "Arm circumference", BMPSITHT = "Sitting height", BMPKNEE = "Knee height", BMPBMI = "Body mass index", PEP12A1 = "Bioelectrical impedance resistance",
  FPP1090 = "Retinal hard exudates", FPP1100 = "Retinal soft exudates", FPP1070 = "Diabetic retinopathy level", FPP1080 = "Retinal hemorrhages or microaneurysms", FPPSURET = "Diabetic retinopathy summary score",
  SPPPEAK = "Peak expiratory flow", RWP = "Red cell distribution width", LMP = "Lymphocyte count", LMPPCNT = "Lymphocyte percentage", LMPDIF = "Lymphocytes (differential)", GRPPCNT = "Granulocyte percentage",
  HAZMNK5R = "Diastolic blood pressure", HAZMNK1R = "Systolic blood pressure", UBP = "Urinary albumin", ATPSI = "Alanine aminotransferase", TIP = "Total iron-binding capacity", PVPSI = "Mean platelet volume",
  AMP = "Serum albumin", CHP = "Serum cholesterol", TCP = "Serum cholesterol (second assay)", AAP = "Apolipoprotein AI", MIP = "Microcytosis", CRP = "C-reactive protein", CEP = "Serum creatinine", SGP = "Serum glucose", GHP = "Glycated haemoglobin",
  HAB4 = "Hospital stays in the past 12 months", HAB5 = "Doctor visits in the past 12 months",
  MAPD_story_recall_delayed_total = "Story recall, delayed", MAPB_story_recall_immediate_total = "Story recall, immediate", HAA_orientation_score = "Orientation (address items)",
  HAT1S = "Walking a mile without stopping, times per month", HAM9S = "Weight at age 25", VEP = "Serum vitamin E", VAP = "Serum vitamin A", LYP = "Serum lycopene", MAPF37 = "Ever used oral estrogen",
  PFP_gait_speed = "Gait speed (8-ft walk)", PFP_steps_8ft_mean = "Steps to walk 8 ft", pulse_mean = "Pulse rate", HAP18_serial_subtraction = "Serial subtraction score",
  HAP17_word_registration = "Three-word registration score", HAP19_word_recall = "Three-word recall score", HAE8D = "Advised to take medication for high cholesterol",
  HAN5AS = "Beans and lentils, times per month", HAN1IS = "Cheese dishes, times per month", HAN7AS = "Margarine, times per month", PBP = "Blood lead")
nice <- function(v, lab) ifelse(v %in% names(NICE), NICE[v], tidy(lab))
