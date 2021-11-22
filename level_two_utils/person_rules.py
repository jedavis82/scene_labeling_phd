from collections import defaultdict
from animal import AnimalRules
from appliances import AppliancesRules
from clothing import ClothingRules
from electronics import ElectronicsRules
from food import FoodRules
from furniture import FurnitureRules
from household import HouseholdRules
from urban import UrbanRules
from vehicle import VehicleRules
from sports import SportsRules


class PersonRules:
    def __init__(self, ontology_df=None, show_sim_result=False):
        """
        Construct the person domain rule bases and compute the person domain level two summary interactions
        :param ontology_df: The object ontology pandas data frame
        :param show_sim_result: If True, show the simulation result for each level two summary. Default False
        """
        if ontology_df is None:
            raise Exception('Must supply ontology data frame')

        self.__ontology_df = ontology_df
        self.__ontology_df.drop(self.__ontology_df[self.__ontology_df['object'] == 'person'].index, inplace=True)
        self.__general_categories_lookup = defaultdict(list)
        self.__domain_categories_lookup = defaultdict(list)
        self.__subdomain_categories_lookup = defaultdict(list)

        self.__create_general_categories_lookup()
        self.__create_domain_categories_lookup()
        self.__create_subdomain_categories_lookup()

        self.__animal_rules = AnimalRules(show_sim_result)
        self.__appliance_rules = AppliancesRules(show_sim_result)
        self.__clothing_rules = ClothingRules(show_sim_result)
        self.__electronics_rules = ElectronicsRules(show_sim_result)
        self.__food_rules = FoodRules(show_sim_result)
        self.__furniture_rules = FurnitureRules(show_sim_result)
        self.__household_rules = HouseholdRules(show_sim_result)
        self.__urban_rules = UrbanRules(show_sim_result)
        self.__vehicle_rules = VehicleRules(show_sim_result)
        self.__sports_rules = SportsRules(show_sim_result)

    def __create_general_categories_lookup(self):
        general_categories = self.__ontology_df['general_category'].unique().flatten().tolist()
        for gc in general_categories:
            objects = list(self.__ontology_df.loc[self.__ontology_df['general_category'] == gc]['object'])
            objects = ['_'.join(o.split(' ')) for o in objects]
            self.__general_categories_lookup[gc] = objects

    def __create_domain_categories_lookup(self):
        domain_categories = self.__ontology_df['subcategory_1'].unique().flatten().tolist()
        for dc in domain_categories:
            objects = list(self.__ontology_df.loc[self.__ontology_df['subcategory_1'] == dc]['object'])
            objects = ['_'.join(o.split(' ')) for o in objects]
            self.__domain_categories_lookup[dc] = objects

    def __create_subdomain_categories_lookup(self):
        subdomain_df = self.__ontology_df.dropna()
        subdomain_categories = subdomain_df['subcategory_2'].unique().flatten().tolist()
        for sc in subdomain_categories:
            objects = list(subdomain_df.loc[subdomain_df['subcategory_2'] == sc]['object'])
            objects = ['_'.join(o.split(' ')) for o in objects]
            self.__subdomain_categories_lookup[sc] = objects

    def __get_general_category(self, label=None):
        for k, v in self.__general_categories_lookup.items():
            if label in v:
                return k

    def __get_domain_category(self, label=None):
        for k, v in self.__domain_categories_lookup.items():
            if label in v:
                return k

    def __get_subdomain_category(self, label=None):
        for k, v in self.__subdomain_categories_lookup.items():
            if label in v:
                return k
        return 'None'  # There is no subdomain for this object

    def get_categories(self, label=None):
        base_label = '_'.join(label.split('_')[:-1])
        general_category = self.__get_general_category(base_label)
        domain_category = self.__get_domain_category(base_label)
        subdomain_category = self.__get_subdomain_category(base_label)
        return general_category, domain_category, subdomain_category, base_label

    def compute_interactions(self, label=None, giou=None, iou=None, sr_angle=None, metadata=None):
        """
        This function will take the input from the various algorithms, compute the general, domain, and subdomain
        categories. Then based on those categories, call the appropriate Fuzzy Control system to determine the human
        and x object interaction.
        metadata will only be supplied if the ref label is "sports ball"
        """
        assert label is not None, 'Must supply a valid ref label'
        assert giou is not None, 'Must supply a valid giou score'
        assert iou is not None, 'Must supply a valid iou score'
        assert sr_angle is not None, 'Must supply a valid spatial relationship angle'
        gen_cat, dom_cat, sub_cat, base_label = self.get_categories(label)

        if gen_cat == 'animal':
            res_label = self.__animal_rules.compute_interaction(base_label, dom_cat, sub_cat, giou, iou, sr_angle)
            return res_label
        elif gen_cat == 'appliances':
            res_label = self.__appliance_rules.compute_interaction(giou, iou, sr_angle)
            return res_label
        elif gen_cat == 'clothing':
            res_label = self.__clothing_rules.compute_interaction(base_label, dom_cat, sub_cat, giou, iou, sr_angle)
            return res_label
        elif gen_cat == 'electronics':
            res_label = self.__electronics_rules.compute_interaction(base_label, dom_cat, sub_cat, giou, iou, sr_angle)
            return res_label
        elif gen_cat == 'food and tableware':
            res_label = self.__food_rules.compute_interaction(base_label, dom_cat, sub_cat, giou, iou, sr_angle)
            return res_label
        elif gen_cat == 'furniture and home decor':
            res_label = self.__furniture_rules.compute_interaction(base_label, dom_cat, sub_cat, giou, iou, sr_angle)
            return res_label
        elif gen_cat == 'household items':
            res_label = self.__household_rules.compute_interaction(base_label, dom_cat, sub_cat, giou, iou, sr_angle)
            return res_label
        elif gen_cat == 'urban':
            res_label = self.__urban_rules.compute_interaction(base_label, dom_cat, sub_cat, giou, iou, sr_angle)
            return res_label
        elif gen_cat == 'vehicle':
            res_label = self.__vehicle_rules.compute_interaction(base_label, dom_cat, sub_cat, giou, iou, sr_angle)
            return res_label
        elif gen_cat == 'sports':
            res_label = self.__sports_rules.compute_interaction(base_label, dom_cat, sub_cat, giou, iou, sr_angle,
                                                                metadata)
            return res_label
        else:
            return None  # Not a valid category. Return None
