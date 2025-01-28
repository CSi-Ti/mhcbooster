

import pandas as pd

from peptdeep.hla.hla_class1 import HLA1_Binding_Classifier
from src.utils.constants import EPSILON
from src.predictors.base_predictor_helper import BasePredictorHelper


class PeptDeepHelper(BasePredictorHelper):
    def __init__(self,
                 peptides: list[str],
                 raw_data: pd.DataFrame,
                 report_directory: str):
        super().__init__('AlphaPeptDeep', report_directory)

        # Prepare dataframe for PeptDeep prediction
        self.peptide_df = pd.DataFrame(peptides, columns=['sequence'])
        self.peptide_df['mods'] = ['' for _ in range(len(self.peptide_df))]
        self.peptide_df['mod_sites'] = ['' for _ in range(len(self.peptide_df))]
        self.peptide_df['charge'] = None
        for col in raw_data.columns:
            if col.startswith('Charge') or col.startswith('charge'):
                self.peptide_df.loc[raw_data[col] == '1', 'charge'] = col[-1]
        self.peptide_df['charge'] = self.peptide_df['charge'].astype(int)
        self.peptide_df['index'] = self.peptide_df.index

    def predict_df(self):
        print('Running AlphaPeptDeep')
        # Predict binding without transfer learning
        model = HLA1_Binding_Classifier()
        model.load_pretrained_hla_model()
        pred_df = model.predict(self.peptide_df)

        # model_mgr = ModelManager()
        # model_mgr.load_installed_models(model_type='HLA')
        #
        # rt_model = model_mgr.rt_model
        # pred_df = rt_model.predict(self.peptide_df)

        # ccs_model = model_mgr.ccs_model
        # ccs_model.predict(pred_df)

        pred_df.sort_values(by=['index'], ascending=True, inplace=True)
        self.pred_df = pred_df
        return self.pred_df

    def score_df(self) -> pd.DataFrame:

        predictions = pd.DataFrame()
        predictions['peptdeep_hla_prob'] = self.pred_df['HLA_prob_pred'].clip(lower=EPSILON).to_numpy()
        # predictions['peptdeep_rt_prob'] = self.pred_df['rt_pred'].clip(lower=EPSILON).to_numpy()
        # predictions['peptdeep_ccs_prob'] = self.pred_df['ccs_pred'].clip(lower=EPSILON).to_numpy()
        # rt_max = np.max(feature_matrix['retentiontime'])
        # rt_delta = np.abs(self.pred_df['rt_pred'] * rt_max - feature_matrix['retentiontime']).clip(lower=EPSILON).to_numpy()
        # predictions['peptdeep_rt_delta'] = rt_delta
        # predictions['peptdeep_rel_rt_delta'] = rt_delta / rt_max

        return predictions

